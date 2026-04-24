/**
 * QSP standalone simulator (model-agnostic driver).
 *
 * Runs the consumer-provided QSP ODE system from a parameter XML and
 * writes the species trajectory to disk. Supports two output formats:
 *
 *   - CSV (human-readable, SimBiology-comparable). Column 0 is Time (days),
 *     remaining columns are species values in original (source) units.
 *   - Raw binary (compact, fast to parse from Python). Used for parameter
 *     sweeps where CSV parse overhead dominates wall time.
 *
 * Binary format (little-endian, packed, no padding):
 *   uint32  magic    = 0x51535042          // "QSPB"
 *   uint32  version  = 1
 *   uint64  n_times                         // number of time rows (incl. t=0)
 *   uint64  n_species                       // columns per row
 *   float64 dt_days
 *   float64 t_end_days
 *   float64 data[n_times * n_species]       // row-major; row = timepoint
 *
 * Time column is NOT stored — it is reconstructible as i*dt for i in
 * [0, n_times), with the final row being t_end_days (possibly clipped).
 *
 * Usage:
 *   qsp_sim --param <xml> [--csv-out <path>] [--binary-out <path>]
 *           [--species-out <path>] [--t-end-days N] [--dt-days N]
 *           [--scenario <scenario.yaml> --drug-metadata <drug_meta.yaml>]
 *
 * With --scenario, doses declared in the scenario YAML are applied as boluses
 * to their target species at exactly their dose times. The solver uses a
 * segmented sampling path: simOdeSample runs between dose boundaries
 * (Nordsieck history preserved, ~5-10x faster than per-step reinit) and
 * setupSamplingRun re-inits CVODE only at dose times. With no doses the
 * segmented loop collapses to a single segment identical to the plain
 * fast-sampling path.
 *
 * Legacy positional form (kept for back-compat with existing tests):
 *   qsp_sim <param_xml> <csv_out> [t_end_days] [dt_days]
 */
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

// Consumer-side generated headers. qsp-codegen emits these into the
// consumer repo; the consumer's CMakeLists adds the directory holding
// them to the include path.
#include "ODE_system.h"
#include "QSPParam.h"
#include "QSP_enum.h"

// Model-init hook. Default no-op ships in default_hooks.cpp inside the
// qsp_sim_core static library; the consumer overrides by defining the
// same symbol in an object file compiled directly into the executable
// (strong definition wins over archive member at link time).
#include "qsp_sim_core/model_hooks.h"

using namespace CancerVCT;

namespace {

struct Args {
    std::string param_file;
    std::string csv_out;
    std::string binary_out;
    std::string species_out;
    std::string compartments_out;  // --compartments-out <path>: one
                                   // compartment name per line, in the
                                   // same order they appear after the
                                   // species block in the v2 binary.
    std::string rules_out;         // --rules-out <path>: one assignment-
                                   // rule name per line, in the same
                                   // order they appear after the
                                   // compartment block in the v2 binary.
    std::string scenario_yaml;
    std::string drug_meta_yaml;
    std::string healthy_yaml;      // --evolve-to-diagnosis <path>: evolve
                                   // from healthy IC; outputs start at t=0
                                   // = diagnosis time. Dose times shift
                                   // by t_diag.
    std::string dump_state_path;   // --dump-state <path>: after evolve,
                                   // serialize ODE state to <path> and exit
                                   // 0. Requires --evolve-to-diagnosis.
    std::string initial_state_path; // --initial-state <path>: skip evolve,
                                   // load ODE state from <path>, set
                                   // t_offset from its header. Mutually
                                   // exclusive with --evolve-to-diagnosis.
    std::string params_hash;       // --params-hash <hex>: stored into
                                   // dumped-state headers; verified against
                                   // --initial-state headers. Catches
                                   // cache/theta mismatches.
    double t_end_days = 365.0;
    double dt_days = 0.1;
};

// Evolve-cache file format (QSTH blob). Fixed 128-byte header followed by
// a CVODEBase full-state payload (see CVODEBase::saveFullState).
//   uint32   magic         = 0x53545148  ('QSTH' little-endian)
//   uint32   version       = 1
//   uint64   n_species_var                (matches current ODE_system on load)
//   float64  t_diagnosis_days             (model-time at which evolve stopped)
//   float64  vt_diameter_cm               (diagnostic V_T diameter at dump)
//   char[32] params_hash_hex              (null-padded; empty means unchecked)
//   char[64] reserved                     (zero-filled for forward compat)
constexpr uint32_t QSTH_MAGIC = 0x53545148u;
constexpr uint32_t QSTH_VERSION = 1u;
constexpr size_t QSTH_HEADER_SIZE = 128;
constexpr size_t QSTH_HASH_LEN = 32;

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " --param <xml> [--csv-out <path>] [--binary-out <path>]\n"
        << "                  [--species-out <path>] [--compartments-out <path>]\n"
        << "                  [--rules-out <path>] [--t-end-days N] [--dt-days N]\n"
        << "                  [--scenario <scenario.yaml> --drug-metadata <drug_meta.yaml>]\n"
        << "                  [--evolve-to-diagnosis <healthy_state.yaml>]\n"
        << "                  [--dump-state <path> | --initial-state <path>]\n"
        << "                  [--params-hash <hex>]\n"
        << "\n"
        << "Evolve cache:\n"
        << "  --dump-state    <path>  with --evolve-to-diagnosis: run evolve, write\n"
        << "                          post-evolve ODE state to <path> (QSTH blob),\n"
        << "                          then exit 0. No scenario sim runs.\n"
        << "  --initial-state <path>  skip evolve; load ODE state from <path> and\n"
        << "                          proceed into the scenario sim. t_offset is\n"
        << "                          read from the blob header.\n"
        << "  --params-hash   <hex>   stored in dumps, verified on loads (catches\n"
        << "                          cache/theta mismatches).\n"
        << "\n"
        << "Legacy: " << prog
        << " <param_xml> <csv_out> [t_end_days] [dt_days]\n";
}

bool parse_args(int argc, char* argv[], Args& out) {
    int i = 1;
    // Optional legacy positional form: <param_xml> <csv_out> [t_end] [dt].
    // Consume positionals that don't look like flags (argv[k][0] != '-'),
    // then fall through to the flag loop so callers can mix positional
    // with flags like --scenario.
    if (argc > i && argv[i][0] != '-') { out.param_file = argv[i++]; }
    if (argc > i && argv[i][0] != '-') { out.csv_out    = argv[i++]; }
    if (argc > i && argv[i][0] != '-') { out.t_end_days = std::stod(argv[i++]); }
    if (argc > i && argv[i][0] != '-') { out.dt_days    = std::stod(argv[i++]); }

    for (; i < argc; ++i) {
        std::string a = argv[i];
        auto need_val = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << name << " requires a value" << std::endl;
                return nullptr;
            }
            return argv[++i];
        };
        if (a == "--param") {
            const char* v = need_val("--param"); if (!v) return false;
            out.param_file = v;
        } else if (a == "--csv-out") {
            const char* v = need_val("--csv-out"); if (!v) return false;
            out.csv_out = v;
        } else if (a == "--binary-out") {
            const char* v = need_val("--binary-out"); if (!v) return false;
            out.binary_out = v;
        } else if (a == "--species-out") {
            const char* v = need_val("--species-out"); if (!v) return false;
            out.species_out = v;
        } else if (a == "--compartments-out") {
            const char* v = need_val("--compartments-out"); if (!v) return false;
            out.compartments_out = v;
        } else if (a == "--rules-out") {
            const char* v = need_val("--rules-out"); if (!v) return false;
            out.rules_out = v;
        } else if (a == "--t-end-days") {
            const char* v = need_val("--t-end-days"); if (!v) return false;
            out.t_end_days = std::stod(v);
        } else if (a == "--dt-days") {
            const char* v = need_val("--dt-days"); if (!v) return false;
            out.dt_days = std::stod(v);
        } else if (a == "--scenario") {
            const char* v = need_val("--scenario"); if (!v) return false;
            out.scenario_yaml = v;
        } else if (a == "--drug-metadata") {
            const char* v = need_val("--drug-metadata"); if (!v) return false;
            out.drug_meta_yaml = v;
        } else if (a == "--evolve-to-diagnosis") {
            const char* v = need_val("--evolve-to-diagnosis"); if (!v) return false;
            out.healthy_yaml = v;
        } else if (a == "--dump-state") {
            const char* v = need_val("--dump-state"); if (!v) return false;
            out.dump_state_path = v;
        } else if (a == "--initial-state") {
            const char* v = need_val("--initial-state"); if (!v) return false;
            out.initial_state_path = v;
        } else if (a == "--params-hash") {
            const char* v = need_val("--params-hash"); if (!v) return false;
            out.params_hash = v;
        } else if (a == "-h" || a == "--help") {
            return false;
        } else {
            std::cerr << "Unknown argument: " << a << std::endl;
            return false;
        }
    }

    if (out.param_file.empty()) {
        std::cerr << "--param is required" << std::endl;
        return false;
    }
    // --dump-state runs evolve only, then exits — trajectory output is
    // meaningless in that mode, so skip the csv/binary requirement.
    if (out.dump_state_path.empty()
        && out.csv_out.empty() && out.binary_out.empty()) {
        std::cerr << "At least one of --csv-out or --binary-out is required"
                  << std::endl;
        return false;
    }
    if (!out.scenario_yaml.empty() && out.drug_meta_yaml.empty()) {
        std::cerr << "--scenario requires --drug-metadata" << std::endl;
        return false;
    }
    if (!out.dump_state_path.empty() && out.healthy_yaml.empty()) {
        std::cerr << "--dump-state requires --evolve-to-diagnosis" << std::endl;
        return false;
    }
    if (!out.initial_state_path.empty() && !out.healthy_yaml.empty()) {
        std::cerr << "--initial-state and --evolve-to-diagnosis are mutually "
                     "exclusive" << std::endl;
        return false;
    }
    if (!out.initial_state_path.empty() && !out.dump_state_path.empty()) {
        std::cerr << "--initial-state and --dump-state are mutually exclusive"
                  << std::endl;
        return false;
    }
    return true;
}

// ----- QSTH (evolve-cache) header I/O --------------------------------
//
// Hand-rolled fixed-size header so the Python side can sanity-check the
// cache without instantiating qsp_sim. Keep in sync with
// qsp_hpc/cpp/evolve_cache.py.

struct QsthHeader {
    uint32_t magic;
    uint32_t version;
    uint64_t n_species_var;
    double t_diagnosis_days;
    double vt_diameter_cm;
    char params_hash[QSTH_HASH_LEN];  // null-padded ASCII hex
    // char reserved[64] — implicit padding; written as zero bytes below.
};

void write_qsth_header(std::ostream& os,
                       uint64_t n_species_var,
                       double t_diagnosis_days,
                       double vt_diameter_cm,
                       const std::string& params_hash) {
    QsthHeader h{};
    h.magic = QSTH_MAGIC;
    h.version = QSTH_VERSION;
    h.n_species_var = n_species_var;
    h.t_diagnosis_days = t_diagnosis_days;
    h.vt_diameter_cm = vt_diameter_cm;
    const size_t copy_n = std::min(params_hash.size(), QSTH_HASH_LEN);
    std::memcpy(h.params_hash, params_hash.data(), copy_n);

    os.write(reinterpret_cast<const char*>(&h.magic), sizeof(h.magic));
    os.write(reinterpret_cast<const char*>(&h.version), sizeof(h.version));
    os.write(reinterpret_cast<const char*>(&h.n_species_var),
             sizeof(h.n_species_var));
    os.write(reinterpret_cast<const char*>(&h.t_diagnosis_days),
             sizeof(h.t_diagnosis_days));
    os.write(reinterpret_cast<const char*>(&h.vt_diameter_cm),
             sizeof(h.vt_diameter_cm));
    os.write(h.params_hash, QSTH_HASH_LEN);
    // Pad to QSTH_HEADER_SIZE.
    const size_t written = sizeof(h.magic) + sizeof(h.version)
        + sizeof(h.n_species_var) + sizeof(h.t_diagnosis_days)
        + sizeof(h.vt_diameter_cm) + QSTH_HASH_LEN;
    static_assert(sizeof(uint32_t) * 2 + sizeof(uint64_t)
                  + sizeof(double) * 2 + QSTH_HASH_LEN <= QSTH_HEADER_SIZE,
                  "QSTH header fields exceed fixed header size");
    const size_t pad = QSTH_HEADER_SIZE - written;
    static const char zeros[QSTH_HEADER_SIZE] = {0};
    os.write(zeros, static_cast<std::streamsize>(pad));
}

// Returns the parsed header. Throws on magic/version mismatch or truncation.
QsthHeader read_qsth_header(std::istream& is) {
    QsthHeader h{};
    is.read(reinterpret_cast<char*>(&h.magic), sizeof(h.magic));
    is.read(reinterpret_cast<char*>(&h.version), sizeof(h.version));
    is.read(reinterpret_cast<char*>(&h.n_species_var),
            sizeof(h.n_species_var));
    is.read(reinterpret_cast<char*>(&h.t_diagnosis_days),
            sizeof(h.t_diagnosis_days));
    is.read(reinterpret_cast<char*>(&h.vt_diameter_cm),
            sizeof(h.vt_diameter_cm));
    is.read(h.params_hash, QSTH_HASH_LEN);
    const size_t read_n = sizeof(h.magic) + sizeof(h.version)
        + sizeof(h.n_species_var) + sizeof(h.t_diagnosis_days)
        + sizeof(h.vt_diameter_cm) + QSTH_HASH_LEN;
    char discard[QSTH_HEADER_SIZE];
    is.read(discard, static_cast<std::streamsize>(QSTH_HEADER_SIZE - read_n));
    if (!is) {
        throw std::runtime_error(
            "QSTH header truncated — expected " + std::to_string(QSTH_HEADER_SIZE)
            + " bytes");
    }
    if (h.magic != QSTH_MAGIC) {
        std::ostringstream msg;
        msg << "QSTH magic mismatch: got 0x" << std::hex << h.magic
            << ", expected 0x" << QSTH_MAGIC
            << " (file is not an evolve-state blob)";
        throw std::runtime_error(msg.str());
    }
    if (h.version != QSTH_VERSION) {
        std::ostringstream msg;
        msg << "QSTH version mismatch: got " << h.version
            << ", expected " << QSTH_VERSION
            << " (rebuild the evolve cache against the current qsp_sim)";
        throw std::runtime_error(msg.str());
    }
    return h;
}

// ---- Species-name → index -------------------------------------------
// getHeader() returns "V_C.nCD4,V_C.Treg,...". The Nth name's index in
// _species_var is N. Only a handful of drug-target species are looked up,
// so a linear scan is fine.
int species_index(const std::string& name) {
    static std::vector<std::string> names;
    if (names.empty()) {
        std::string h = ODE_system::getHeader();
        std::string cur;
        for (char c : h) {
            if (c == ',') { names.push_back(cur); cur.clear(); }
            else          { cur += c; }
        }
        if (!cur.empty()) names.push_back(cur);
    }
    for (size_t i = 0; i < names.size(); i++) {
        if (names[i] == name) return static_cast<int>(i);
    }
    return -1;
}

// ---- Dose scheduling ------------------------------------------------
struct Bolus {
    double t_sim;        // absolute simulation time in solver units
    int species_idx;
    double amount;       // in the storage unit of the species (SI moles, cells, etc.)
    std::string label;   // for logging
};

double scale_amount_to_storage(const std::string& units, double amount_in_units) {
    // Storage is SI substance. The wrapper applies doses as "add N mol to the
    // amount" rather than "set concentration", mirroring MATLAB sbiodose with
    // AmountUnits='mole'.
    if (units == "mole")      return amount_in_units;
    if (units == "cell")      return amount_in_units / 6.02214076e23;
    // Some species (e.g. cyclophosphamide's V_C.Cy) are stored natively
    // in milligrams in the SBML (substanceUnits=milligram); an mg dose
    // is added to the state directly.
    if (units == "milligram") return amount_in_units;
    throw std::runtime_error("unknown dose units: " + units);
}

// Expand a scenario + drug metadata into a list of Bolus events.
// time_factor converts days to the solver's time unit.
std::vector<Bolus> build_dose_plan(
    const YAML::Node& scenario,
    const YAML::Node& drug_meta,
    double time_factor)
{
    std::vector<Bolus> plan;
    if (!scenario["dosing"]) return plan;
    const auto& dosing = scenario["dosing"];

    const double patient_weight = dosing["patientWeight"]
        ? dosing["patientWeight"].as<double>() : 70.0;
    const double patient_bsa = dosing["patientBSA"]
        ? dosing["patientBSA"].as<double>() : 1.9;

    if (!dosing["drugs"]) return plan;
    for (const auto& drug_n : dosing["drugs"]) {
        std::string drug = drug_n.as<std::string>();

        const auto& meta_drugs = drug_meta["drugs"];
        if (!meta_drugs || !meta_drugs[drug]) {
            throw std::runtime_error("drug not in drug_metadata.yaml: " + drug);
        }
        const auto& md = meta_drugs[drug];
        const std::string units = md["units"].as<std::string>();
        const std::string basis = md["dose_basis"].as<std::string>();

        const std::string dose_key = drug + "_dose";
        const std::string sched_key = drug + "_schedule";
        if (!dosing[dose_key]) {
            throw std::runtime_error("scenario is missing " + dose_key);
        }
        const double raw_dose = dosing[dose_key].as<double>();
        const auto sched = dosing[sched_key].as<std::vector<double>>();
        if (sched.size() != 3) {
            throw std::runtime_error(sched_key + " must be [start, interval, repeat]");
        }
        const double start_day = sched[0];
        const double interval_day = sched[1];
        const int repeat = static_cast<int>(sched[2]);

        double total_amount = 0.0;
        if (basis == "per_weight") {
            const double mw = md["mw"].as<double>();
            total_amount = patient_weight * raw_dose / mw;
        } else if (basis == "per_bsa") {
            const double mw = md["mw"].as<double>();
            total_amount = patient_bsa * raw_dose / mw;
        } else if (basis == "direct") {
            total_amount = raw_dose;
        } else {
            throw std::runtime_error("unknown dose_basis: " + basis);
        }

        for (const auto& target : md["targets"]) {
            const std::string sp_name = target["species"].as<std::string>();
            const double frac = target["fraction"].as<double>();
            const int idx = species_index(sp_name);
            if (idx < 0) {
                throw std::runtime_error("target species not found in ODE: " + sp_name);
            }
            const double storage_amount = scale_amount_to_storage(
                units, total_amount * frac);

            for (int r = 0; r < repeat; r++) {
                const double t_day = start_day + r * interval_day;
                plan.push_back({
                    t_day * time_factor,
                    idx,
                    storage_amount,
                    drug + "@" + sp_name,
                });
            }
        }
    }
    return plan;
}

}  // namespace

int main(int argc, char* argv[]) {
    Args args;
    if (!parse_args(argc, argv, args)) {
        print_usage(argv[0]);
        return 1;
    }

#ifdef MODEL_UNITS
    const double time_factor = 1.0;
#else
    const double time_factor = 86400.0;
#endif
    double t_end = args.t_end_days * time_factor;
    const double dt = args.dt_days * time_factor;

    QSPParam param;
    param.initializeParams(args.param_file);
    ODE_system::setup_class_parameters(param);

    ODE_system ode;
    ode.setup_instance_variables(param);
    ode.setup_instance_tolerance(param);
    ode.eval_init_assignment();

    const std::string header = ODE_system::getHeader();
    const size_t n_species = static_cast<size_t>(ode.getNumOutputSpecies());

    if (!args.species_out.empty()) {
        std::ofstream sp_out(args.species_out);
        size_t start = 0;
        for (size_t i = 0; i <= header.size(); ++i) {
            if (i == header.size() || header[i] == ',') {
                sp_out << header.substr(start, i - start) << '\n';
                start = i + 1;
            }
        }
    }

    // Compartment volumes + non-compartment assignment-rule values are
    // emitted alongside species so downstream calibration-target functions
    // can read derived quantities (V_T, phi_collagen, C_total, …) by name
    // without re-deriving them in Python. MATLAB SimBiology already exposes
    // these in `simdata.Data`; matching the column set keeps the C++
    // backend a drop-in replacement. Names come from the codegen so they
    // stay in sync with the SBML.
    const std::vector<std::string> extra_comps = ODE_system::getCompartmentNames();
    const std::vector<std::string> extra_rules = ODE_system::getAssignmentRuleNames();
    const size_t n_compartments = extra_comps.size();
    const size_t n_rules = extra_rules.size();

    if (!args.compartments_out.empty()) {
        std::ofstream c_out(args.compartments_out);
        for (const auto& c : extra_comps) c_out << c << '\n';
    }
    if (!args.rules_out.empty()) {
        std::ofstream r_out(args.rules_out);
        for (const auto& r : extra_rules) r_out << r << '\n';
    }

    std::ofstream csv;
    if (!args.csv_out.empty()) {
        csv.open(args.csv_out);
        csv << std::scientific << std::setprecision(12);
        // operator<< on CVODEBase emits all species prefixed with commas; pair
        // it with "Time,<header>" (no leading comma in getHeader()) for a
        // self-consistent CSV.
        csv << "Time," << header;
        for (const auto& c : extra_comps) csv << "," << c;
        for (const auto& r : extra_rules) csv << "," << r;
        csv << std::endl;
    }

    // Binary v2 layout (header is 56 bytes):
    //   uint32 magic, uint32 version=2, uint64 n_t,
    //   uint64 n_species, uint64 n_compartments, uint64 n_rules,
    //   double dt_days, double t_end_days
    // followed by n_t × (n_sp + n_comp + n_rules) doubles in row-major
    // order (species first, then compartments, then rules — matching the
    // order of the *_out name files). v1 (40-byte header, species-only
    // body) is no longer produced; readers should accept it in
    // backward-compat mode but can't be generated here anymore.
    //
    // Header is written with a placeholder n_t; we seek back and patch
    // it after stepping because float dt may yield one more or one fewer
    // step than ceil(t_end/dt).
    std::ofstream bin;
    const uint32_t MAGIC = 0x51535042u;  // "QSPB"
    const uint32_t VERSION = 2;
    if (!args.binary_out.empty()) {
        bin.open(args.binary_out, std::ios::binary);
        uint64_t n_times_placeholder = 0;
        uint64_t n_sp64 = static_cast<uint64_t>(n_species);
        uint64_t n_comp64 = static_cast<uint64_t>(n_compartments);
        uint64_t n_rules64 = static_cast<uint64_t>(n_rules);
        bin.write(reinterpret_cast<const char*>(&MAGIC), sizeof(MAGIC));
        bin.write(reinterpret_cast<const char*>(&VERSION), sizeof(VERSION));
        bin.write(reinterpret_cast<const char*>(&n_times_placeholder), sizeof(uint64_t));
        bin.write(reinterpret_cast<const char*>(&n_sp64), sizeof(uint64_t));
        bin.write(reinterpret_cast<const char*>(&n_comp64), sizeof(uint64_t));
        bin.write(reinterpret_cast<const char*>(&n_rules64), sizeof(uint64_t));
        bin.write(reinterpret_cast<const char*>(&args.dt_days), sizeof(double));
        bin.write(reinterpret_cast<const char*>(&args.t_end_days), sizeof(double));
    }

    const size_t n_cols = n_species + n_compartments + n_rules;
    std::vector<double> row(n_cols);

    // Optional: replace ICs with the healthy microinvasive state and integrate
    // forward until V_T diameter crosses the target. The state at return is
    // the "diagnosis" state; CSV output starts at user-time 0 (i.e. diagnosis).
    // t_offset is the solver's internal time at diagnosis. The scenario
    // dosing loop runs from t_offset to t_offset + t_end while write_state
    // writes (t - t_offset) so the output time axis stays user-relative.
    //
    // Three entry paths:
    //   (a) --evolve-to-diagnosis alone: integrate, then continue into the
    //       scenario sim.
    //   (b) --evolve-to-diagnosis + --dump-state: integrate, write QSTH blob,
    //       exit 0. The same theta can then be run under many scenarios via
    //       (c) without re-paying the evolve cost.
    //   (c) --initial-state: skip evolve entirely, load state from a prior
    //       dump, set t_offset from the blob header. Scenario runs as usual.
    double t_offset = 0.0;
    if (!args.healthy_yaml.empty()) {
        EvolveOpts eo;
        eo.yaml_path = args.healthy_yaml;
        eo.time_factor = time_factor;
        eo.verbose = true;
        EvolveResult er = evolve_to_diagnosis(ode, eo);
        if (!er.success) {
            std::cerr << "evolve_to_diagnosis REJECTED: "
                      << er.reject_reason << std::endl;
            return 2;
        }
        std::cerr << "evolve_to_diagnosis: t_diag=" << er.t_diagnosis_days
                  << " d, diameter=" << er.diameter_cm << " cm\n";
        t_offset = er.t_diagnosis_days * time_factor;

        if (!args.dump_state_path.empty()) {
            std::ofstream dump(args.dump_state_path, std::ios::binary);
            if (!dump) {
                std::cerr << "failed to open --dump-state path: "
                          << args.dump_state_path << std::endl;
                return 3;
            }
            write_qsth_header(dump,
                              static_cast<uint64_t>(ode.get_num_variables()),
                              er.t_diagnosis_days,
                              er.diameter_cm,
                              args.params_hash);
            ode.saveFullState(dump);
            if (!dump) {
                std::cerr << "failed to write --dump-state payload: "
                          << args.dump_state_path << std::endl;
                return 3;
            }
            dump.close();
            std::cerr << "dump-state: wrote post-evolve ODE state to "
                      << args.dump_state_path << " (t_diag="
                      << er.t_diagnosis_days << " d)\n";
            return 0;
        }
    } else if (!args.initial_state_path.empty()) {
        std::ifstream in(args.initial_state_path, std::ios::binary);
        if (!in) {
            std::cerr << "failed to open --initial-state path: "
                      << args.initial_state_path << std::endl;
            return 3;
        }
        QsthHeader h{};
        try {
            h = read_qsth_header(in);
        } catch (const std::exception& e) {
            std::cerr << "initial-state header error: " << e.what() << "\n";
            return 3;
        }
        if (!args.params_hash.empty()) {
            const std::string stored(h.params_hash,
                strnlen(h.params_hash, QSTH_HASH_LEN));
            if (stored != args.params_hash) {
                std::cerr << "initial-state params_hash mismatch: file='"
                          << stored << "' arg='" << args.params_hash
                          << "' (cache and current theta don't match — "
                             "rebuild the cache)" << std::endl;
                return 3;
            }
        }
        try {
            ode.loadFullState(in);
        } catch (const std::exception& e) {
            std::cerr << "initial-state payload error: " << e.what() << "\n";
            return 3;
        }
        // Sync CVODE's N_Vector with the freshly-loaded _species_var and
        // refresh _species_other. The subsequent setupSamplingRun at
        // t=t_offset will CVodeReInit on top of this, so CVODE's internal
        // step history is discarded (correct — we have no history from the
        // evolve run we skipped).
        ode.updateVar();
        t_offset = h.t_diagnosis_days * time_factor;
        std::cerr << "initial-state: loaded ODE state from "
                  << args.initial_state_path
                  << " (t_diag=" << h.t_diagnosis_days
                  << " d, diameter=" << h.vt_diameter_cm << " cm)\n";
    }

    auto write_state = [&](double t) {
        if (csv.is_open()) {
            csv << (t - t_offset) / time_factor << ode;
            for (const auto& c : extra_comps) {
                csv << "," << ode.get_compartment_volume(c);
            }
            for (const auto& r : extra_rules) {
                csv << "," << ode.get_assignment_rule_value(r);
            }
            csv << std::endl;
        }
        if (bin.is_open()) {
            for (size_t i = 0; i < n_species; ++i) {
                row[i] = ode.getSpeciesOutputValue(static_cast<int>(i));
            }
            for (size_t i = 0; i < n_compartments; ++i) {
                row[n_species + i] = ode.get_compartment_volume(extra_comps[i]);
            }
            for (size_t i = 0; i < n_rules; ++i) {
                row[n_species + n_compartments + i] =
                    ode.get_assignment_rule_value(extra_rules[i]);
            }
            bin.write(reinterpret_cast<const char*>(row.data()),
                      static_cast<std::streamsize>(n_cols * sizeof(double)));
        }
    };

    // Load scenario + drug metadata if requested, building the bolus plan
    // before any integration so we know whether to take the fast sampling
    // path or the step-based dosing path.
    std::vector<Bolus> dose_plan;
    if (!args.scenario_yaml.empty()) {
        YAML::Node scenario = YAML::LoadFile(args.scenario_yaml);
        YAML::Node drug_meta = YAML::LoadFile(args.drug_meta_yaml);
        dose_plan = build_dose_plan(scenario, drug_meta, time_factor);
        // Shift into solver time: user-specified dose times are relative to
        // diagnosis (t=0 user = t_offset solver).
        for (auto& b : dose_plan) b.t_sim += t_offset;
        std::cerr << "Loaded " << dose_plan.size() << " bolus events from "
                  << args.scenario_yaml << std::endl;
        for (const auto& b : dose_plan) {
            std::cerr << "  t=" << (b.t_sim - t_offset) / time_factor << "d  "
                      << b.label << "  amount=" << b.amount << std::endl;
        }
        if (scenario["sim_config"] && scenario["sim_config"]["stop_time"]) {
            t_end = scenario["sim_config"]["stop_time"].as<double>() * time_factor;
        }
    }

    const double t_stop = t_offset + t_end;

    // Segmented sampling: partition [t_offset, t_stop] on dose boundaries
    // and use simOdeSample within each segment, so CVODE's Nordsieck history
    // is preserved across output ticks. Only re-init (setupSamplingRun)
    // happens at dose times. With no doses, this collapses to a single
    // segment and matches the previous fast-sampling path exactly.
    //
    // Also gives doses exact timing: a bolus at day 7.2 is applied at 7.2
    // even if output ticks land at 7.0 and 7.5 — the old step-based path
    // integrated through dose times and applied at the next tick boundary,
    // which was only accidentally correct when doses aligned with ticks.
    std::vector<double> dose_boundaries;
    {
        std::vector<double> tmp;
        for (const auto& b : dose_plan) {
            if (b.t_sim > t_offset && b.t_sim <= t_stop) tmp.push_back(b.t_sim);
        }
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        dose_boundaries = std::move(tmp);
    }

    auto apply_at = [&](double dose_t) {
        for (const auto& b : dose_plan) {
            if (b.t_sim != dose_t) continue;
            double cur = ode.getSpeciesVar(
                static_cast<unsigned int>(b.species_idx), false);
            ode.setSpeciesVar(
                static_cast<unsigned int>(b.species_idx),
                cur + b.amount, false);
            std::cerr << "  [dose] t=" << (dose_t - t_offset) / time_factor
                      << "d  " << b.label << "  +" << b.amount << std::endl;
        }
    };

    // Apply any t=0 (solver = t_offset) boluses before integrating so the
    // first written sample reflects post-dose state.
    apply_at(t_offset);

    uint64_t n_times = 1;
    write_state(t_offset);

    double t = t_offset;
    size_t next_dose_idx = 0;
    // Tick index counts emitted output points (0 = the t=t_offset row above).
    // Recomputing next_tick = t_offset + tick_idx * dt each iteration avoids
    // accumulated float drift over many steps — important because we compare
    // it to seg_end with exact equality.
    int tick_idx = 1;

    double seg_end = dose_boundaries.empty()
        ? t_stop : std::min(dose_boundaries.front(), t_stop);
    ode.setupSamplingRun(seg_end, t);

    int step = 0;
    while (t < t_stop) {
        double next_tick = t_offset + tick_idx * dt;
        if (next_tick > t_stop) next_tick = t_stop;
        double target = std::min(next_tick, seg_end);

        ode.simOdeSample(target);
        t = target;

        // Exact float comparisons are safe here: `target` was assigned from
        // one of {next_tick, seg_end, t_stop} and simOdeSample returns the
        // caller's t verbatim. next_tick is recomputed from tick_idx*dt, not
        // accumulated.
        bool at_boundary = (t == seg_end) && (next_dose_idx < dose_boundaries.size());
        bool at_tick = (t == next_tick);

        // Apply bolus before the output emit so the written sample is
        // post-dose (matching the step-based M10 port's semantics and the
        // MATLAB sbiosimulate dose_schedule semantics).
        if (at_boundary) apply_at(t);

        if (at_tick) {
            write_state(t);
            n_times++;
            tick_idx++;
        }

        if (at_boundary && t < t_stop) {
            ++next_dose_idx;
            seg_end = (next_dose_idx < dose_boundaries.size())
                ? std::min(dose_boundaries[next_dose_idx], t_stop) : t_stop;
            ode.setupSamplingRun(seg_end, t);
        }

        ++step;
        if (step % 1000 == 0) {
            std::cerr << "  t=" << (t - t_offset) / time_factor
                      << " days" << std::endl;
        }
    }

    if (csv.is_open()) {
        csv.close();
        std::cerr << "Wrote " << n_times << " time points to " << args.csv_out
                  << std::endl;
    }
    if (bin.is_open()) {
        // n_times slot is at offset 8 in the header (after magic + version).
        bin.seekp(sizeof(uint32_t) * 2, std::ios::beg);
        bin.write(reinterpret_cast<const char*>(&n_times), sizeof(uint64_t));
        bin.close();
        std::cerr << "Wrote " << n_times << " time points × " << n_cols
                  << " columns ("
                  << n_species << " species + "
                  << n_compartments << " compartments + "
                  << n_rules << " rules) to "
                  << args.binary_out << std::endl;
    }
    return 0;
}