% bench_matlab_cold.m - Single load+simulate, used by the wall-clock probe.
%
% Variables expected:
%   sbml_path  : path to SBML
%   stop_time  : integration end time (days)
%
% Time is captured by the *parent* Python script around the matlab -batch
% invocation; this script just performs one full load+simulate so that
% wall-clock measurement reflects what an SBI workflow pays per cold call.

required = {'sbml_path', 'stop_time'};
for i = 1:numel(required)
    if ~exist(required{i}, 'var')
        error('bench_matlab_cold:missingVar', '%s not set', required{i});
    end
end

model = sbmlimport(sbml_path);
cfg = getconfigset(model);
set(cfg, 'SolverType', 'sundials');
cfg.SolverOptions.RelativeTolerance = 1e-6;
cfg.SolverOptions.AbsoluteTolerance = 1e-9;
cfg.StopTime = stop_time;
cfg.SolverOptions.OutputTimes = (0:0.1:stop_time)';

sbiosimulate(model);
