% build_model.m - Build a multi-compartment SimBiology model and export it.
%
% Run once to (re)generate model.sbml. Both the MATLAB benchmark loop and
% the qsp-codegen pipeline read this file, so both sides solve the same
% problem.
%
% Usage:
%   matlab -batch "run('paper/benchmark/build_model.m')"
%
% Model: a chain of N compartments connected by bidirectional first-order
% transfers with parameters spanning two orders of magnitude (multi-
% timescale stiffness), plus a Michaelis-Menten clearance term in every
% compartment. An initial bolus lives in compartment 1.
%
% Why this shape: a small toy PK model is too easy — both engines bottom
% out on output-writing rather than integration. With N≈25 species and a
% sparse, near-tridiagonal Jacobian, the analytical-Jacobian + KLU sparse
% solver path that qsp-codegen takes diverges meaningfully from
% SimBiology's default finite-difference + dense LU, while remaining
% small enough to ship in this directory and rebuild quickly.

here = fileparts(mfilename('fullpath'));
out_sbml = fullfile(here, 'model.sbml');

N = 25;
rng(42);  % reproducible parameter draws

model = sbiomodel('chain_pbpk');

% Compartments (modest volume variation so initial concentrations differ).
comps = cell(N, 1);
for i = 1:N
    vol = 1.0 + 0.5 * rand();
    comps{i} = addcompartment(model, sprintf('c%d', i), vol);
    addspecies(comps{i}, 'Drug', 0.0);
end

% Bolus in compartment 1.
comps{1}.Species(1).InitialAmount = 100.0;

% Bidirectional transfer between adjacent compartments. Rates drawn over
% [0.1, 10] /day so the system has both fast equilibria and slow
% drift -- this is what makes CVODE work.
for i = 1:(N - 1)
    kf = 10 ^ (2 * rand() - 1);  % 0.1 .. 10
    kb = 10 ^ (2 * rand() - 1);
    addparameter(model, sprintf('kf_%d', i), kf);
    addparameter(model, sprintf('kb_%d', i), kb);

    rf = addreaction(model, sprintf('c%d.Drug -> c%d.Drug', i, i + 1));
    klf = addkineticlaw(rf, 'MassAction');
    set(klf, 'ParameterVariableNames', {sprintf('kf_%d', i)});

    rb = addreaction(model, sprintf('c%d.Drug -> c%d.Drug', i + 1, i));
    klb = addkineticlaw(rb, 'MassAction');
    set(klb, 'ParameterVariableNames', {sprintf('kb_%d', i)});
end

% Michaelis-Menten clearance per compartment. Vmax / Km drawn so a few
% compartments saturate quickly while others stay near-linear -- sympy
% emits a non-trivial diagonal of Jacobian entries that depend on the
% local Drug level.
for i = 1:N
    vmax = 10 ^ (3 * rand() - 2);  % 0.01 .. 10
    km   = 10 ^ (2 * rand() - 1);  % 0.1  .. 10
    addparameter(model, sprintf('Vmax_%d', i), vmax);
    addparameter(model, sprintf('Km_%d', i),   km);

    rcl = addreaction(model, sprintf('c%d.Drug -> null', i));
    set(rcl, 'ReactionRate', sprintf( ...
        'Vmax_%d * c%d.Drug / (Km_%d + c%d.Drug)', i, i, i, i));
end

cfg = getconfigset(model);
set(cfg, 'SolverType', 'sundials');
cfg.SolverOptions.RelativeTolerance = 1e-6;
cfg.SolverOptions.AbsoluteTolerance = 1e-9;
cfg.StopTime = 365;
cfg.SolverOptions.OutputTimes = (0:0.1:365)';

sbmlexport(model, out_sbml);
fprintf('Wrote %s (N=%d compartments, %d reactions, %d parameters)\n', ...
        out_sbml, N, length(model.Reactions), length(model.Parameters));
