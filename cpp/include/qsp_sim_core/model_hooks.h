#ifndef QSP_SIM_CORE_MODEL_HOOKS_H
#define QSP_SIM_CORE_MODEL_HOOKS_H

// Model-init hook interface the generic qsp_sim driver calls into.
//
// The driver does not know how a specific model (PDAC, HCC, ...) sets up
// its "initial" state. Consumers who need a pre-scenario evolve phase
// (e.g. evolve from healthy tissue until a tumor diameter is reached)
// provide an implementation of `evolve_to_diagnosis` and link it into
// the executable. Consumers that don't need this behavior pick up the
// default no-op shipped in `default_hooks.cpp` (link order wins via the
// standard "strong definition overrides archive member" rule; the
// default sits in a static library so any consumer-provided object file
// with the same symbol supplants it).
//
// The hook is a free function rather than a virtual on ODE_system to
// keep the generated-code ABI untouched: qsp-codegen emits ODE_system,
// and layering a virtual there would require touching the generator on
// every model-init design change.

#include <string>
#include <vector>

namespace CancerVCT {

class ODE_system;

struct EvolveOpts {
    std::string yaml_path;          // Model-specific evolve spec. Semantics
                                    // are the consumer's responsibility;
                                    // the driver only forwards the path.
    double target_diameter_cm = -1; // Tumor-model convention; unused by
                                    // models that don't track a diameter.
    double tumor_cells = -1;
    double time_factor = 86400.0;   // 86400 (SI sec) or 1.0 (model-unit days)
    bool verbose = false;

    // Optional: dump dense state during the burn-in phase. When
    // ``trajectory_path`` is non-empty, the consumer's
    // ``evolve_to_diagnosis`` implementation should open a binary v2
    // file at this path (same magic 0x51535042 / VERSION 2 layout the
    // post-scenario writer in qsp_sim_main produces) and append a row
    // every ``trajectory_dt_days`` of model time. ``extra_comps`` /
    // ``extra_rules`` carry the column layout (compartment volumes
    // and assignment-rule values to track in addition to species),
    // matching whatever the post-scenario CLI is configured to write
    // so downstream readers can use one schema. Use 0.0 for
    // ``trajectory_dt_days`` to fall back to the spec's step_days.
    std::string trajectory_path;
    double trajectory_dt_days = 0.0;
    std::vector<std::string> trajectory_extra_comps;
    std::vector<std::string> trajectory_extra_rules;
};

struct EvolveResult {
    bool success = false;
    double t_diagnosis_days = 0.0;  // Returned model-time (days) at which
                                    // the consumer considers the state
                                    // "ready" for the scenario sim. The
                                    // driver shifts dose times by this.
    double diameter_cm = 0.0;
    std::string reject_reason;
};

// Consumer override point. Default (in default_hooks.cpp) returns
// success=true, t_diagnosis_days=0 — i.e. no pre-scenario evolve.
EvolveResult evolve_to_diagnosis(ODE_system& ode, const EvolveOpts& opts);

}  // namespace CancerVCT

#endif
