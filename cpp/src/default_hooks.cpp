// Default (no-op) implementation of the model-init hooks declared in
// qsp_sim_core/model_hooks.h. Consumers override by linking their own
// definition into the executable: since this default lives in the
// qsp_sim_core static library, an object file in the consumer target
// that defines `evolve_to_diagnosis` takes precedence at link time.
//
// A no-op evolve leaves the ODE at its post-eval_init_assignment state
// and reports t_diagnosis_days = 0, which is the right behavior for
// models whose scenario sims start directly from ICs.

#include "qsp_sim_core/model_hooks.h"

namespace qsp_sim_core {

EvolveResult evolve_to_diagnosis(ODE_system& /*ode*/, const EvolveOpts& /*opts*/) {
    EvolveResult r;
    r.success = true;
    r.t_diagnosis_days = 0.0;
    r.diameter_cm = 0.0;
    return r;
}

}  // namespace qsp_sim_core
