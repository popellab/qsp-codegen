#ifndef QSP_SIM_CORE_TRAJECTORY_WRITER_H
#define QSP_SIM_CORE_TRAJECTORY_WRITER_H

// v3 trajectory binary writer (the "QSPB" format).
//
// Owns the byte-level layout so consumers — both qsp_sim_main and
// out-of-tree drivers like pdac-build's evolve_to_diagnosis — write the
// same bytes through one code path. When the format bumps (v3→v4),
// editing this header is the whole change; consumers that don't rebuild
// against the new qsp_sim_core get a clean compile/link mismatch instead
// of silently producing unreadable binaries (the v2→v3 drift this helper
// was extracted to prevent — see GH#18).
//
// Layout (little-endian, packed, no padding; 80-byte header):
//   uint32  magic              = 0x51535042   // "QSPB"
//   uint32  version            = 3
//   uint64  n_times                            // patched by finalize()
//   uint64  n_species
//   uint64  n_compartments
//   uint64  n_rules
//   float64 min_cadence_hours
//   float64 t_end_days
//   float64 t_offset_days                      // patched by finalize()
//   uint64  n_cvode_steps                      // patched by finalize()
//   uint64  reserved           = 0
//   float64 data[n_times * (1 + n_species + n_compartments + n_rules)]
//
// Each row is [t_user_days, species..., compartments..., rules...].
// Sample times are non-uniform: CV_ONE_STEP cadence with a floor (see
// qsp_sim_main.cpp's main loop), which is why v3 stores time per-row.
//
// Usage:
//   std::ofstream os(path, std::ios::binary);
//   TrajectoryWriter w(os, n_species, n_compartments, n_rules,
//                      min_cadence_hours, t_end_days);
//   for (each sample) w.write_row(t_days, row_ptr, n_state_cols);
//   w.finalize(t_offset_days, n_cvode_steps);
//
// The stream MUST be opened in binary mode and MUST be seekable
// (finalize() rewinds to patch header fields).

#include <cstdint>
#include <cstring>
#include <ostream>
#include <stdexcept>

namespace CancerVCT {

class TrajectoryWriter {
public:
    static constexpr uint32_t MAGIC = 0x51535042u;  // "QSPB"
    static constexpr uint32_t VERSION = 3u;
    static constexpr size_t HEADER_SIZE = 80;

    TrajectoryWriter(std::ostream& os,
                     uint64_t n_species,
                     uint64_t n_compartments,
                     uint64_t n_rules,
                     double min_cadence_hours,
                     double t_end_days)
        : os_(os),
          n_species_(n_species),
          n_compartments_(n_compartments),
          n_rules_(n_rules),
          n_state_cols_(n_species + n_compartments + n_rules)
    {
        const uint32_t magic = MAGIC;
        const uint32_t version = VERSION;
        const uint64_t n_times_placeholder = 0;
        const double t_offset_days_placeholder = 0.0;
        const uint64_t n_cvode_steps_placeholder = 0;
        const uint64_t reserved = 0;

        write_raw(&magic, sizeof(magic));
        write_raw(&version, sizeof(version));
        write_raw(&n_times_placeholder, sizeof(uint64_t));
        write_raw(&n_species_, sizeof(uint64_t));
        write_raw(&n_compartments_, sizeof(uint64_t));
        write_raw(&n_rules_, sizeof(uint64_t));
        write_raw(&min_cadence_hours, sizeof(double));
        write_raw(&t_end_days, sizeof(double));
        write_raw(&t_offset_days_placeholder, sizeof(double));
        write_raw(&n_cvode_steps_placeholder, sizeof(uint64_t));
        write_raw(&reserved, sizeof(uint64_t));
    }

    // Append one row. `state_cols` is [species..., compartments..., rules...]
    // with length n_species + n_compartments + n_rules; `n_state_cols`
    // is passed in for caller-side defensive checking.
    void write_row(double t_days, const double* state_cols, size_t n_state_cols) {
        if (n_state_cols != n_state_cols_) {
            throw std::runtime_error(
                "TrajectoryWriter::write_row: state column count mismatch");
        }
        write_raw(&t_days, sizeof(double));
        write_raw(state_cols, n_state_cols * sizeof(double));
        ++n_times_;
    }

    // Patch the three deferred header fields (n_times, t_offset_days,
    // n_cvode_steps). Idempotent: safe to call once at end-of-run.
    void finalize(double t_offset_days, uint64_t n_cvode_steps) {
        constexpr std::streamoff OFF_N_TIMES = 8;        // after magic+version
        constexpr std::streamoff OFF_T_OFFSET = 56;      // after t_end_days
        constexpr std::streamoff OFF_N_CVODE_STEPS = 64; // after t_offset_days
        os_.seekp(OFF_N_TIMES, std::ios::beg);
        os_.write(reinterpret_cast<const char*>(&n_times_), sizeof(uint64_t));
        os_.seekp(OFF_T_OFFSET, std::ios::beg);
        os_.write(reinterpret_cast<const char*>(&t_offset_days), sizeof(double));
        os_.seekp(OFF_N_CVODE_STEPS, std::ios::beg);
        os_.write(reinterpret_cast<const char*>(&n_cvode_steps), sizeof(uint64_t));
    }

    uint64_t n_times() const { return n_times_; }
    size_t n_state_cols() const { return n_state_cols_; }

private:
    void write_raw(const void* p, size_t n) {
        os_.write(reinterpret_cast<const char*>(p),
                  static_cast<std::streamsize>(n));
    }

    std::ostream& os_;
    uint64_t n_species_;
    uint64_t n_compartments_;
    uint64_t n_rules_;
    size_t n_state_cols_;
    uint64_t n_times_ = 0;
};

}  // namespace CancerVCT

#endif  // QSP_SIM_CORE_TRAJECTORY_WRITER_H
