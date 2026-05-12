#include "qsp_sim_core/CVODEBase.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

const int mxstep = 500000;

CVODEBase::CVODEBase()
: _species_var()
, _species_other()
, _nonspecies_var()
, _neq(0)
, _y(NULL)
, _nroot(0)
, _nevent(0)
, _delayEvents()
, _A(NULL)
, _LS(NULL)
, _cvode_mem(NULL)
, _sunctx(NULL)
, _trigger_element_type()
, _trigger_element_satisfied()
, _event_triggered()
{
	SUNContext_Create(SUN_COMM_NULL, &_sunctx);
}

CVODEBase::CVODEBase(const CVODEBase & c)
: _species_var(c._species_var)
, _species_other(c._species_other)
, _nonspecies_var(c._nonspecies_var)
, _neq(c._neq)
, _y(NULL)
, _nroot(0)
, _nevent(0)
, _delayEvents()
, _A(NULL)
, _LS(NULL)
, _cvode_mem(NULL)
, _sunctx(NULL)
, _trigger_element_type()
, _trigger_element_satisfied()
, _event_triggered()
{
	SUNContext_Create(SUN_COMM_NULL, &_sunctx);
}
CVODEBase::~CVODEBase()
{
	freeMem();
	if (_sunctx) { SUNContext_Free(&_sunctx); _sunctx = NULL; }
}

/*! simulate model for some time
\param [in] tStart: simulation starting time
\param [in] tStep: simulation interval

Before every sub step simulated, the solver is reset
due to the need of updating species concentration from
PDE module. While doing this, t and species are also reset
to accommodate the need to deserialize.

Similarly, after the interval, results are copied to vector
version of container for easier serialization and manipulation
by other modules.

TODO:
The error CV_TOO_CLOSE ("tout too close to t0 to start integration")
is issued if one of the following occurs:
(1) tout == t0, or
(2) |tout - t0| < 2 * eps * max(|t0|, |tout|)
In the future this function may need to be updated so that t0 is set
to 0 every time step, with solver reconstructed instead of reinitiated.
*/
void CVODEBase::simOdeStep(double tStart, double tStep){

	int flag;

	realtype t = tStart;
	realtype tEnd = tStart + tStep;
	realtype t1 = tEnd;

	//std::cout << "beginning: " << t << ", " << tEnd << std::endl;

	// in case events need to be executed at the beginning
	// of a step.
	// This could happen when:
	// 1. t = 0
	// 2. ODE modified externally between steps.
	resolveEvents(t);

	while (t < tEnd){

		bool delayedExecution = false;
		bool discontinuity = false;
		t1 = tEnd;

		realtype tNextDisc = 0;
		bool queueNotEmpty = getNexTimeDisc(tNextDisc);
		//std::cout << "delayed execution time: " << tNextDisc 
		//	<< ", delay: " << queueNotEmpty << std::endl;
		if (queueNotEmpty && tNextDisc < tEnd)
		{
			t1 = tNextDisc;
			delayedExecution = true;
		}
		//std::cout << "cycle: " << t << ", " << t1 << std::endl;

		resetSolver(t, t1);

		flag = CVode(_cvode_mem, t1, _y, &t, CV_NORMAL);

		if (flag == CV_TOO_CLOSE)
		{	// just move forward in this case
			t = t1;
		}
		else {
			check_flag(&flag, "CVode", 1);
		}

		if (flag == CV_ROOT_RETURN)
		{
			int* rootsFound = new int[_nroot];
			CVodeGetRootInfo(_cvode_mem, rootsFound);

			// update variable value triggered events
			updateTriggerComponentConditionsOnRoot(rootsFound);

			discontinuity = evaluateAllEvents(t);
			delete[] rootsFound;
		}
		else {
			/**/
			if (delayedExecution) {
				// update trigger components for non-persistent events
				// this is actually not necessary. we keep it there 
				// just in case.
				updateTriggerComponentConditionsOnValue(t);

				//std::cout << "delay queue size: " << _delayEvents.size() << std::endl;
				int eventToExecute = _delayEvents.back().second;
				realtype ttemp;
				bool delay = eventExecution(eventToExecute, true, ttemp);
				_delayEvents.pop_back();
				//std::cout << "delay queue size: " << _delayEvents.size() << std::endl;
				discontinuity = !delay;
			}
		}

		// either delayed execution or root found 
		if (discontinuity)
		{
			resolveEvents(t);
			resetTransient();
			/*
			std::cout << "t_off: " << _nonspecies_var[3] << std::endl;

			std::cout << "Trigger: ";
			for (auto i = 0; i < _nroot; i++)
			{
				std::cout << _trigger_element_satisfied[i] << ",";
			}
			std::cout << std::endl;*/
		}

		save_y();
		update_y_other();
	}
	std::cout << std::flush;
	std::cerr << std::flush;
	//std::cout << "End of step:" << t << *this << std::endl;
}

void CVODEBase::resolveEvents(realtype t){
	bool discontinuity = true;
	while (discontinuity){
		updateTriggerComponentConditionsOnValue(t);
		discontinuity = evaluateAllEvents(t);
	}
}

/*! Manually update variable values.
	Need to do this after altering variable values
	externally and need to get output before simulating 
	any time step.
*/
void CVODEBase::updateVar(void)
{
	restore_y();
	update_y_other();
}

/*! Setup the solver
This Should be called only once during construction, when no
prior allocation of memory block to solver has taken place.

In this function, serial type containers are constructed;
variable values are copied to serial type container;
memory block is allocated to solver, which is initiated to t=0;
Other settings are configured such as tolerance, user data
and linear solver.

*/
void CVODEBase::setupCVODE(){

	bool res = true;

	_neq = _species_var.size();

	try{
		int flag;

		_y = N_VNew_Serial(_neq, _sunctx);
		check_flag((void *)_y, "N_VNew_Serial", 0);

		//_abstol = N_VNew_Serial(_neq, _sunctx);
		//check_flag((void *)_abstol, "N_VNew_Serial", 0);

		_cvode_mem = CVodeCreate(CV_BDF, _sunctx);
		check_flag((void *)_cvode_mem, "CVodeCreate", 0);

		/* Call CVodeInit to initialize the integrator memory and specify the
		* user's right hand side function in y'=f(t,y), the inital time T0, and
		* the initial dependent variable vector y.
		* need to do this in derived class: f and g are defined in derived
		* class; they are static functions, which cannot be virtual functions
		* so we cannot declare them in the base class.
		*/

		initSolver(0);


		flag = CVodeSetMaxNumSteps(_cvode_mem, mxstep);

		/* Enforce non-negativity: all species are biological quantities
		 * (concentrations, cell counts, synapse complexes) that cannot
		 * be negative.  Without this, CVODE can transiently drive species
		 * negative during Newton iteration, producing NaN in expressions
		 * like pow(phi_collagen, 1.8) and collapsing h to machine zero. */
		{
			N_Vector constraints = N_VNew_Serial(_neq, _sunctx);
			for (int i = 0; i < _neq; i++) {
				NV_DATA_S(constraints)[i] = 1.0;  // y[i] >= 0
			}
			flag = CVodeSetConstraints(_cvode_mem, constraints);
			N_VDestroy(constraints);
		}

		/* Passing the pointer of this system to the solver, so that
		* the parameters etc. can be accessed from the static function f and g.
		* Function f () and g (root finding) need to access object
		* specific version of parameter values.
		* this has some cost on performance. */
		flag = CVodeSetUserData(_cvode_mem, this);
		check_flag(&flag, "CVodeSVtolerances", 1);

		/* Linear solver selection: if the derived class exposes an
		 * analytical sparse Jacobian (nnz > 0) AND KLU is compiled in,
		 * use SUNSparseMatrix + SUNLinSol_KLU. Otherwise fall back to
		 * the dense path. Sparse is ~2-5x faster per integration step on
		 * the PDAC model (164 species, ~5% density) since it skips the
		 * O(n^3) dense LU in favour of KLU's symbolic + numeric
		 * factorization over the static sparsity pattern.
		 */
		sunindextype jac_nnz = getJacobianNnz();
		CVLsJacFn jac_fn = getJacobianFn();
#ifdef USE_KLU
		if (jac_nnz > 0 && jac_fn != nullptr) {
			_A = SUNSparseMatrix(_neq, _neq, jac_nnz, CSC_MAT, _sunctx);
			check_flag(&flag, "SUNSparseMatrix", 1);
			_LS = SUNLinSol_KLU(_y, _A, _sunctx);
			check_flag(&flag, "SUNLinSol_KLU", 1);
			flag = CVodeSetLinearSolver(_cvode_mem, _LS, _A);
			check_flag(&flag, "CVodeSetLinearSolver (KLU)", 1);
			flag = CVodeSetJacFn(_cvode_mem, jac_fn);
			check_flag(&flag, "CVodeSetJacFn", 1);
		} else
#endif
		{
			/* Dense fallback (always available). We do NOT attach the
			 * codegen's analytical Jacobian here: it is emitted against
			 * SUNSparseMatrix (CSC), not SUNDenseMatrix, so using it on
			 * the dense path would segfault in SUNSparseMatrix_IndexPointers.
			 * CVODE's default finite-difference Jacobian is correct here.
			 */
			_A = SUNDenseMatrix(_neq, _neq, _sunctx);
			check_flag(&flag, "SUNDenseMatrix", 1);
			_LS = SUNLinSol_Dense(_y, _A, _sunctx);
			check_flag(&flag, "SUNLinSol_Dense", 1);
			flag = CVodeSetLinearSolver(_cvode_mem, _LS, _A);
			check_flag(&flag, "CVodeSetLinearSolver (dense)", 1);
		}

	}
	catch (std::string s){
		std::cerr << "Initiating CVODE solver, error: " << s;
		exit(1);
	}
	//std::cout << "_neq: " << _neq << std::endl;
}

/*! get the t of the next potential discontinuity in simulation.
This can be:
1. Delay of execution from variable-associated trigger condition
(detected with root finding function)
Time-associated trigger condition in one event used to be handled here.
Now they are delt with in rootfinding function.
*/
bool CVODEBase::getNexTimeDisc(realtype& t){
	/**/
	if (!_delayEvents.empty())
	{
		t = _delayEvents.back().first;
		return true;
	}
	else {
		return false;
	}
	return false;
}
/* update trigger conditions
*/
void CVODEBase::updateTriggerComponentConditionsOnRoot(int* rootsFound) {
	for (auto i = 0; i < _nroot; i++)
	{
		if (isTransient(i))
		{
			if (rootsFound[i] != 0)
			{
				_trigger_element_satisfied[i] = isTransientEq(i);
			}
		}
		else {
			if (rootsFound[i] == 1 && !_trigger_element_satisfied[i])
			{
				_trigger_element_satisfied[i] = true;
			}
			else if (rootsFound[i] == -1 && _trigger_element_satisfied[i])
			{
				_trigger_element_satisfied[i] = false;
			}
		}
	}
	return;
}

void CVODEBase::updateTriggerComponentConditionsOnValue(realtype t) {
	for (auto i = 0; i < _nroot; i++)
	{
		_trigger_element_satisfied[i] = triggerComponentEvaluate(i, t,
			_trigger_element_satisfied[i]);
	}
	//std::cout << std::endl;
	return;
}

/* evaluate event triggers
*/
bool CVODEBase::evaluateAllEvents(realtype t) {
	// evaluate all events, check if they are up for execution.
	//std::cout << "trigger evaluations: " << std::endl;
	bool exec = false;
	for (auto i = 0; i < _nevent; i++)
	{
		bool trigger = eventEvaluate(i);
		//std::cout << "event " << i <<
		//	", evaluation result: " << trigger << std::endl;
		// only execute when evaluation result change from false to true
		if (trigger && !_event_triggered[i])
		{
			realtype dt = 0;
			bool setDelay = eventExecution(i, false, dt);
			/**/
			if (setDelay)
			{
				_delayEvents.push_back(std::make_pair(t + dt, i));
				std::sort(_delayEvents.rbegin(), _delayEvents.rend());
			}
			//std::cout << "execution: event " << i << std::endl;
			// any event executed?
			exec |= !setDelay;
		}
		_event_triggered[i] = trigger;
	}
	return exec;
}
/*  reset triggers
*/
void CVODEBase::resetEventTriggers() {
	resetTransient();
	for (auto i = 0; i < _nevent; i++)
	{
		_event_triggered[i] = eventEvaluate(i);
	}
	return;
}
/*  reset all transient conditions to false
this is called after event evaluation.
*/
void CVODEBase::resetTransient() {
	for (auto i = 0; i < _nroot; i++)
	{
		if (isTransientEq(i))
		{
			_trigger_element_satisfied[i] = false;
		}
		else {
			if (isTransientNeq(i)) {
				_trigger_element_satisfied[i] = true;
			}
		}
	}
}
/*!
 * Get and print some final statistics
 */
void CVODEBase::PrintFinalStats(void *cvode_mem)
{
	long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
	int flag;

	flag = CVodeGetNumSteps(cvode_mem, &nst);
	check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	check_flag(&flag, "CVodeGetNumRhsEvals", 1);
	flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
	check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	check_flag(&flag, "CVodeGetNumErrTestFails", 1);
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

	flag = CVodeGetNumJacEvals(cvode_mem, &nje);
	check_flag(&flag, "CVDlsGetNumJacEvals", 1);
	flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
	check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

	flag = CVodeGetNumGEvals(cvode_mem, &nge);
	check_flag(&flag, "CVodeGetNumGEvals", 1);

	printf("\nFinal Statistics:\n");
	printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
		nst, nfe, nsetups, nfeLS, nje);
	printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
		nni, ncfn, netf, nge);
}

double CVODEBase::getSpeciesVar(unsigned int idx, bool raw) const
{
	if (idx < _species_var.size()){
		if (raw)
		{
			return _species_var[idx] / get_unit_conversion_species(idx);
		}
		else{
			return _species_var[idx];
		}
	}
	else{
		throw std::invalid_argument("Accessing ODE variable: out of range");
	}
}

void CVODEBase::setSpeciesVar(unsigned int idx, double val, bool raw)
{
	if (idx < _species_var.size()){
		if (raw)
		{
			_species_var[idx] = val * get_unit_conversion_species(idx);
		}
		else{
			_species_var[idx] = val;
		}
	}
	else{
		throw std::invalid_argument("Assignment to ODE variable: out of range");
	}
	return;
}

double CVODEBase::getParameterVal(unsigned int idx, bool raw) const
{
	if (idx < _nonspecies_var.size()){
		if (raw)
		{
			return _nonspecies_var[idx] / get_unit_conversion_nspvar(idx);
		}
		else{
			return _nonspecies_var[idx];
		}
	}
	else{
		throw std::invalid_argument("Accessing ODE instance parameter: out of range");
	}
}

void CVODEBase::setParameterVal(unsigned int idx, double val, bool raw)
{
	if (idx < _nonspecies_var.size()){
		if (raw)
		{
			_nonspecies_var[idx] = val * get_unit_conversion_nspvar(idx);
		}
		else{
			_nonspecies_var[idx] = val;
		}
	}
	else{
		throw std::invalid_argument("Assignment to ODE instance parameter: out of range");
	}
	return;
}

/*! copy _species_var to the serial container _y;
reset start time and initial condition.
\param [in] t0: start time
\param [in] t1: end time
*/
void CVODEBase::resetSolver(realtype t0, realtype t1){
	int flag = 0;
	restore_y();
	flag = CVodeSetStopTime(_cvode_mem, t1);
	flag = CVodeReInit(_cvode_mem, t0, _y);
	return;
}

/*! One-time setup for a sampling run. Syncs _species_var -> _y, resolves
	any initial-assignment events at t=0, pushes that post-init state into
	CVODE via a single CVodeReInit, and pins the stop time so CV_NORMAL
	calls won't walk past tEndOfSim. Must be called once before the first
	simOdeSample call.

	The CVodeReInit here is essential: setupCVODE()/initSolver() was called
	during construction with whatever _y held at that point (typically
	pre-initial-assignment defaults), and CVODE internally caches that
	state. Without the re-init, simOdeSample would integrate from the
	stale state, not the properly-initialized one.

*/
void CVODEBase::setupSamplingRun(double tEndOfSim, double t0, double h0_hint){
	restore_y();
	resolveEvents(t0);
	CVodeReInit(_cvode_mem, t0, _y);
	CVodeSetStopTime(_cvode_mem, tEndOfSim);
	if (h0_hint > 0.0) {
		// Skip CVHin's heuristic h₀ pick. CVHin estimates h₀ from
		// ||f(t0,y0)||_WRMS, which collapses to sub-ULP at solver times
		// O(1e8) s when a bolus jolts a fast-binding subsystem — the
		// step then advances t by less than one double-precision ULP
		// and the outer loop's "did not advance" guard aborts the run.
		// Passing the pre-event h_cur as the hint mirrors SimBiology.
		CVodeSetInitStep(_cvode_mem, static_cast<realtype>(h0_hint));
	}
}

/*! Advance CVODE from its current internal time to tEnd in CV_NORMAL mode.
	Unlike simOdeStep, this does NOT call CVodeReInit — the Nordsieck history
	and adaptive step-size controller are preserved across sampling points,
	which is where the factor-of-5+ speedup over repeated simOdeStep calls
	comes from on stiff systems. CV_TOO_CLOSE is silently treated as a no-op
	(we're already at tEnd).
*/
void CVODEBase::simOdeSample(double tEnd){
	realtype t_ret = tEnd;
	int flag = CVode(_cvode_mem, tEnd, _y, &t_ret, CV_NORMAL);
	if (flag != CV_TOO_CLOSE) {
		check_flag(&flag, "CVode", 1);
	}
	// Sync derived outputs (assignment-rule species, original-unit views)
	// so getSpeciesOutputValue/operator<< emit current state.
	save_y();
	update_y_other();
}

/*! Advance CVODE by one internal step (CV_ONE_STEP) and return the resulting
	integration time. Used by output loops that prefer solver-native cadence
	to a fixed output grid: the caller dumps a row whenever t advances past a
	configurable minimum cadence and at the simulation stop time. CVODE's
	stop time is pinned by the prior setupSamplingRun call, so the step won't
	walk past tEndOfSim regardless of `tEndClamp` (which CV_ONE_STEP uses only
	to determine integration direction).
*/
double CVODEBase::simOdeStepOne(double tEndClamp){
	realtype t_ret = tEndClamp;
	int flag = CVode(_cvode_mem, tEndClamp, _y, &t_ret, CV_ONE_STEP);
	if (flag != CV_TOO_CLOSE) {
		check_flag(&flag, "CVode", 1);
	}
	save_y();
	update_y_other();
	return static_cast<double>(t_ret);
}

/*! Wrapper over CVodeGetNumSteps. Returns 0 on failure. The counter is reset
	by CVodeReInit (so by setupSamplingRun and resetSolver), making it usable
	for per-segment step instrumentation if the caller queries before each
	re-init. Accumulating across segments is the caller's responsibility.
*/
long CVODEBase::getNumSteps() const {
	long n_steps = 0;
	int flag = CVodeGetNumSteps(_cvode_mem, &n_steps);
	if (flag != CV_SUCCESS) {
		return 0;
	}
	return n_steps;
}

CVODEBase::StepStats CVODEBase::getStepStats() const {
	StepStats s{};
	CVodeGetNumSteps(_cvode_mem, &s.nst);
	CVodeGetNumRhsEvals(_cvode_mem, &s.nfe);
	CVodeGetNumErrTestFails(_cvode_mem, &s.netf);
	CVodeGetNumNonlinSolvConvFails(_cvode_mem, &s.ncfn);
	CVodeGetNumJacEvals(_cvode_mem, &s.nje);
	CVodeGetNumNonlinSolvIters(_cvode_mem, &s.nni);
	CVodeGetNumLinSolvSetups(_cvode_mem, &s.nsetups);
	CVodeGetLastOrder(_cvode_mem, &s.last_order);
	realtype hlast = 0.0, hcur = 0.0;
	CVodeGetLastStep(_cvode_mem, &hlast);
	CVodeGetCurrentStep(_cvode_mem, &hcur);
	s.last_h = static_cast<double>(hlast);
	s.cur_h = static_cast<double>(hcur);
	return s;
}
// ----- Evolve-cache full-state serialization -----------------------------
//
// Pure data dump: ODE_system state vectors + delay-event queue + trigger
// flags, nothing CVODE-internal. Used by qsp_sim --dump-state /
// --initial-state so multi-scenario sweeps share one evolve run per theta.

namespace {

template <typename T>
inline void write_pod(std::ostream& os, const T& v) {
    os.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

template <typename T>
inline void read_pod(std::istream& is, T& v) {
    is.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (!is) {
        throw std::runtime_error(
            "CVODEBase::loadFullState: truncated stream");
    }
}

}  // namespace

void CVODEBase::saveFullState(std::ostream& os) const {
    const uint64_t n_sp = static_cast<uint64_t>(_species_var.size());
    const uint64_t n_nsp = static_cast<uint64_t>(_nonspecies_var.size());
    const uint64_t n_de = static_cast<uint64_t>(_delayEvents.size());
    const uint64_t n_tsat = static_cast<uint64_t>(
        _trigger_element_satisfied.size());
    const uint64_t n_etrig = static_cast<uint64_t>(_event_triggered.size());

    write_pod(os, n_sp);
    if (n_sp) {
        os.write(reinterpret_cast<const char*>(_species_var.data()),
                 static_cast<std::streamsize>(n_sp * sizeof(double)));
    }

    write_pod(os, n_nsp);
    if (n_nsp) {
        os.write(reinterpret_cast<const char*>(_nonspecies_var.data()),
                 static_cast<std::streamsize>(n_nsp * sizeof(double)));
    }

    write_pod(os, n_de);
    for (const auto& ev : _delayEvents) {
        const double t = static_cast<double>(ev.first);
        const int32_t idx = static_cast<int32_t>(ev.second);
        write_pod(os, t);
        write_pod(os, idx);
    }

    write_pod(os, n_tsat);
    for (bool b : _trigger_element_satisfied) {
        const uint8_t byte = b ? 1u : 0u;
        write_pod(os, byte);
    }

    write_pod(os, n_etrig);
    for (bool b : _event_triggered) {
        const uint8_t byte = b ? 1u : 0u;
        write_pod(os, byte);
    }
}

void CVODEBase::loadFullState(std::istream& is) {
    uint64_t n_sp = 0, n_nsp = 0, n_de = 0, n_tsat = 0, n_etrig = 0;

    read_pod(is, n_sp);
    if (n_sp != _species_var.size()) {
        std::ostringstream msg;
        msg << "CVODEBase::loadFullState: species_var length mismatch "
            << "(file=" << n_sp << ", current model=" << _species_var.size()
            << "). Likely a qsp_sim version / SBML codegen mismatch — "
               "rebuild the cache.";
        throw std::runtime_error(msg.str());
    }
    if (n_sp) {
        is.read(reinterpret_cast<char*>(_species_var.data()),
                static_cast<std::streamsize>(n_sp * sizeof(double)));
        if (!is) throw std::runtime_error(
            "CVODEBase::loadFullState: truncated species_var payload");
    }

    read_pod(is, n_nsp);
    if (n_nsp != _nonspecies_var.size()) {
        std::ostringstream msg;
        msg << "CVODEBase::loadFullState: nonspecies_var length mismatch "
            << "(file=" << n_nsp
            << ", current model=" << _nonspecies_var.size() << ")";
        throw std::runtime_error(msg.str());
    }
    if (n_nsp) {
        is.read(reinterpret_cast<char*>(_nonspecies_var.data()),
                static_cast<std::streamsize>(n_nsp * sizeof(double)));
        if (!is) throw std::runtime_error(
            "CVODEBase::loadFullState: truncated nonspecies_var payload");
    }

    read_pod(is, n_de);
    _delayEvents.clear();
    _delayEvents.reserve(static_cast<size_t>(n_de));
    for (uint64_t i = 0; i < n_de; ++i) {
        double t = 0.0;
        int32_t idx = 0;
        read_pod(is, t);
        read_pod(is, idx);
        _delayEvents.emplace_back(static_cast<realtype>(t),
                                  static_cast<int>(idx));
    }

    read_pod(is, n_tsat);
    if (n_tsat != _trigger_element_satisfied.size()) {
        std::ostringstream msg;
        msg << "CVODEBase::loadFullState: trigger_element_satisfied length "
               "mismatch (file=" << n_tsat
            << ", current model=" << _trigger_element_satisfied.size() << ")";
        throw std::runtime_error(msg.str());
    }
    for (uint64_t i = 0; i < n_tsat; ++i) {
        uint8_t byte = 0;
        read_pod(is, byte);
        _trigger_element_satisfied[i] = (byte != 0);
    }

    read_pod(is, n_etrig);
    if (n_etrig != _event_triggered.size()) {
        std::ostringstream msg;
        msg << "CVODEBase::loadFullState: event_triggered length mismatch "
            << "(file=" << n_etrig
            << ", current model=" << _event_triggered.size() << ")";
        throw std::runtime_error(msg.str());
    }
    for (uint64_t i = 0; i < n_etrig; ++i) {
        uint8_t byte = 0;
        read_pod(is, byte);
        _event_triggered[i] = (byte != 0);
    }
}

/*! copy variable value from vector to serial
*/
void CVODEBase::restore_y(){
	for (auto i = 0; i < _neq; i++)
	{
		NV_DATA_S(_y)[i] = _species_var[i];
	}
	/* Clear _species_var to save memory
	_species_var.clear(); */
}

/*! copy variable value from serial to vector
*/
void CVODEBase::save_y(){

	//_species_var.resize(_neq, 0);
	for (auto i = 0; i < _neq; i++)
	{
		_species_var[i] = NV_DATA_S(_y)[i] *
			(NV_DATA_S(_y)[i] < 0 ? allow_negative(i) : 1 );
	}
}

/*! get variable value with original unit
*/
realtype CVODEBase::getVarOriginalUnit(int i) const{
	realtype v = 0;
	if (i < _neq){
		v = _species_var[i];
	}
	else{
		v = _species_other[i - _neq];
	}
	v /= get_unit_conversion_species(i);
	return v;
}

/*! free memory blocks.
This include solver memory block
and any other serial type blocks
*/
bool CVODEBase::freeMem(){

	/* Free y and abstol vectors */
	N_VDestroy(_y);
	//N_VDestroy(_abstol);

	/* Free integrator memory */
	CVodeFree(&_cvode_mem);
	return true;
}
/*!
* Check function return value...
*   - opt == 0 means SUNDIALS function allocates memory so check if
*            returned NULL pointer
*   - opt == 1 means SUNDIALS function returns a flag so check if
*            flag >= 0
*   - opt == 2 means function allocates memory so check if returned
*            NULL pointer
*/
void CVODEBase::check_flag(void *flagvalue, const char *funcname, int opt)
{
	int *errflag;
	std::stringstream ss;
	bool res = false;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		ss << "\nSUNDIALS_ERROR: " << funcname << "() failed - returned NULL pointer\n";
		throw ss.str();
	}

	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *)flagvalue;
		if (*errflag < 0) {
			ss << "\nSUNDIALS_ERROR: " << funcname << "() failed with flag = " << *errflag << "\n";
			throw ss.str();
		}
	}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL) {

		ss << "\nMEMORY_ERROR: " << funcname << "() failed - returned NULL pointer\n";
		throw ss.str();
	}
}