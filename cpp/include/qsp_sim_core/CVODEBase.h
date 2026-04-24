#ifndef __CVODE_BASE__
#define __CVODE_BASE__

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/utility.hpp> /* std::pair */

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <cvode/cvode_ls.h>            /* CVLsJacFn type for analytical J hook */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix (KLU)    */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#ifdef USE_KLU
#include <sunlinsol/sunlinsol_klu.h>   /* sparse direct via SuiteSparse KLU    */
#endif
#include <sundials/sundials_types.h>            /* defs. of sunrealtype, sunindextype */
#include <sundials/sundials_types_deprecated.h> /* realtype compat (removed in v7)    */
#include <sundials/sundials_context.h>          /* SUNContext (SUNDIALS >= 6)         */



#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

typedef std::vector< double > state_type;

//! Base class for CVode 
class CVODEBase
{
protected:
	enum EVENT_TRIGGER_ELEM_TYPE {
		TRIGGER_NON_INSTANT,
		TRIGGER_EQ,
		TRIGGER_NEQ
	};

public:
	CVODEBase(); 
	CVODEBase(const CVODEBase & c);
	~CVODEBase();

	//! simulate ODE model (event-handling path, used by ABM coupling).
	//! Calls CVodeReInit at every step boundary — cheap on short steps but
	//! prohibitive when the caller wants fine-grained output sampling of
	//! a long simulation.
	void simOdeStep(double tStart, double tStep);

	//! Forward-advance CVODE to absolute time tEnd without reinitializing.
	//! Designed for output sampling loops in event-free simulations: the
	//! caller sets up once via setupSamplingRun(tEnd_of_simulation) and
	//! then invokes simOdeSample(t_out_i) for each snapshot time. CVODE's
	//! internal adaptive step / Nordsieck history is preserved across
	//! calls, which can be 5-10x faster than simOdeStep on stiff systems
	//! with many output points per event. Does NOT handle SBML events or
	//! delayed-execution discontinuities — use simOdeStep for those.
	void simOdeSample(double tEnd);

	//! One-time setup before the first simOdeSample call: sets the CVODE
	//! stop time to t_end_of_sim so CV_NORMAL stepping won't walk past it,
	//! and runs any pending initial-assignment events at t0.
	//!
	//! `t0` is the absolute integration start time. Pass the ODE's current
	//! internal time when entering a sampling run after a prior
	//! simOdeStep-based phase (e.g. evolve_to_diagnosis); if the fast path
	//! runs from fresh ICs, leave it at the 0.0 default.
	void setupSamplingRun(double tEndOfSim, double t0 = 0.0);

	//! examples of optional output
	void PrintFinalStats(void *cvode_mem);

	//! ODE state to stream
	friend std::ostream & operator<<(std::ostream &os, const CVODEBase & ode) ;

	//! species varaible value with original units
	double getSpeciesVar(unsigned int idx, bool raw = true)const;
	//! set species varaible value with original units
	void setSpeciesVar(unsigned int idx, double val, bool raw = true);
	//! non species varaible value with original units
	double getParameterVal(unsigned int idx, bool raw = true)const;
	//! set non species varaible value with original units
	void setParameterVal(unsigned int idx, double val, bool raw = true);

	//! number of species emitted by operator<< (sp_var + sp_other)
	int getNumOutputSpecies() const {
		return _neq + static_cast<int>(_species_other.size());
	}
	//! value of output species i in original units — public forwarder over
	//! the protected virtual getVarOriginalUnit, for callers that need to
	//! snapshot the full row to a buffer rather than to an ostream
	double getSpeciesOutputValue(int i) const { return getVarOriginalUnit(i); }

	//! Analytical Jacobian hooks for sparse-KLU linear solver.
	//! Default implementation returns 0 / nullptr, which keeps setupCVODE
	//! on the dense linear-solver path (backwards compatible). A derived
	//! class with a generated ODE_system::jac overrides these to expose
	//! the sparse Jacobian; setupCVODE then allocates a SUNSparseMatrix
	//! of the right nnz and registers the callback via CVodeSetJacFn.
	virtual sunindextype getJacobianNnz() const { return 0; }
	virtual CVLsJacFn getJacobianFn() const { return nullptr; }

	//! manually update solver variable values
	void updateVar(void);

	//! Serialize the full instance-level ODE state (species + nonspecies +
	//! delay-event queue + trigger/event flags) to a binary stream. Used
	//! by the evolve-to-diagnosis cache (qsp_sim --dump-state) to preserve
	//! the post-evolve state across scenario runs. Does NOT capture CVODE
	//! internal state (Nordsieck history etc.) — after loadFullState the
	//! caller must sync via updateVar() and re-init the solver via
	//! setupSamplingRun() or resetSolver() before integrating.
	//!
	//! Field-by-field layout (little-endian, packed), read by loadFullState:
	//!   uint64 n_species_var;      float64[n_species_var]
	//!   uint64 n_nonspecies_var;   float64[n_nonspecies_var]
	//!   uint64 n_delay_events;     { float64 t; int32 idx; }[n_delay_events]
	//!   uint64 n_trigger_sat;      uint8[n_trigger_sat]    (0/1)
	//!   uint64 n_event_trig;       uint8[n_event_trig]     (0/1)
	void saveFullState(std::ostream& os) const;
	//! Mirror of saveFullState: read the same layout and populate the
	//! instance state vectors in place. Throws std::runtime_error on
	//! length/stream errors. Caller is responsible for calling updateVar()
	//! + setupSamplingRun() / resetSolver() afterwards.
	void loadFullState(std::istream& is);


protected:

	//! initial setup, prepare memory blocks for the solver
	void setupCVODE();
	//! pure virtual. Pass rhs function and initial conditions to solver, instantiated in derived class
	virtual void initSolver(realtype t0) = 0;
	//! setup variables
	virtual void setupVariables(void) = 0;
	//! setup events 
	virtual void setupEvents(void) = 0;
	//! reset starting time and _y after model is already initialted
	void resetSolver(realtype t0, realtype t1);
	//! copy _species_var and save to _y
	void restore_y();
	//! copy _y and save to _species_var
	void save_y();
	//! update species that are not parf of lhs of ode
	virtual void update_y_other(void) = 0;
	//! evaluate one trigger component 
    virtual bool triggerComponentEvaluate(int i, realtype t, bool curr) = 0;
	//! update trigger conditions when root found
	void updateTriggerComponentConditionsOnRoot(int* rootsFound);
	//! update trigger conditions at time t, without root 
	void updateTriggerComponentConditionsOnValue(realtype t);
	//! evaluate event triggers
	bool evaluateAllEvents(realtype t);
	//! resolve event assignments recursively
	void resolveEvents(realtype t);
	//! reset trigger after evaluation
	void resetEventTriggers();
	//! reset transient trigger conditions to false
	void resetTransient();
	//! trigger condition associated with one root 
	bool getSatisfied(int i);
	//! if trigger condition is '==' or '!='
	bool isTransient(int i);
	//! if trigger condition is '==' 
	bool isTransientEq(int i);
	//! if trigger condition is '!='
	bool isTransientNeq(int i);
	//! evaluate one event trigger
	virtual bool eventEvaluate(int i) = 0;
	//! execute one event
	virtual bool eventExecution(int i, bool delay, realtype& dt) = 0;
	//! get variable value with original unit
	virtual double getVarOriginalUnit(int i) const;
	//! get unit conversion scalor
	virtual realtype get_unit_conversion_species(int i) const = 0;
	//! get unit conversion scalor
	virtual realtype get_unit_conversion_nspvar(int i) const = 0;
	//! check if a variable is allowed to become negative
	virtual bool allow_negative(int i) const {return true;};


	//! some functions defined by SBML interpretor
	inline static double root(double a, double b) { return std::pow(b, 1.0 / a); };

	//! check the return values of CVode related functions
	void check_flag(void *flagvalue, const char *funcname, int opt);

	//! variable species. Species in the left-hand side of ODEs 
	state_type _species_var;
	//! other speceis. Listed as species in SBML, but not in lhs. 
	state_type _species_other;
	//! Non-species variable subject to event assignments 
	state_type _nonspecies_var;
	//! constant species/parameters/compartments 
	//state_type _parameter_const;

	//! number of equations (same as nr of _species_var)
	int _neq;
	//! pointer to variable species, in format compatible with CVode
	N_Vector _y;
	//! number of rootfinding functions
	int _nroot;
	//! number of events
	int _nevent;

	//! delayed events sorted vector
	std::vector<std::pair <realtype, int> > _delayEvents;

	//! SUNMatrix for linear solver 
	SUNMatrix _A;
	//! linear solver object for CVode
	SUNLinearSolver _LS;
	//! solver memory block
	void * _cvode_mem;
	//! SUNDIALS context (required by SUNDIALS >= 6)
	SUNContext _sunctx;

	//! event only triggered when g(y, t) = 0. Serialization not needed. 
	std::vector<EVENT_TRIGGER_ELEM_TYPE>  _trigger_element_type;
	//! one event trigger element is satisfied 
	std::vector<bool>  _trigger_element_satisfied;
	//! one event is triggered
	std::vector<bool>  _event_triggered;


private:
	friend class boost::serialization::access;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);

	bool freeMem();
	//! get next t for potential discontinuity.
	bool getNexTimeDisc(realtype& t);



};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(CVODEBase)

inline std::ostream & operator<<(std::ostream &os, const CVODEBase & ode){
	int nrSpecies = ode._neq + ode._species_other.size();
	for (auto i = 0; i < nrSpecies; i++)
	{
		os << "," << ode.getVarOriginalUnit(i);
	}
	return os;
}

inline bool CVODEBase::getSatisfied(int i) {
	return _trigger_element_satisfied[i];
}

inline bool CVODEBase::isTransient(int i) {
	return (_trigger_element_type[i] == TRIGGER_NON_INSTANT ? false: true);
}
inline bool CVODEBase::isTransientEq(int i) {
	return (_trigger_element_type[i] == TRIGGER_EQ ? true: false);
}
inline bool CVODEBase::isTransientNeq(int i) {
	return (_trigger_element_type[i] == TRIGGER_NEQ ? true: false);
}
template<class Archive>
inline void CVODEBase::serialize(Archive & ar, const unsigned int /* version */){
	ar & BOOST_SERIALIZATION_NVP(_species_var);
	ar & BOOST_SERIALIZATION_NVP(_nonspecies_var);
	ar & BOOST_SERIALIZATION_NVP(_delayEvents);
	ar & BOOST_SERIALIZATION_NVP(_trigger_element_satisfied);
	ar & BOOST_SERIALIZATION_NVP(_event_triggered);
}

#endif
