% bench_matlab.m - Time SimBiology integration for the benchmark harness.
%
% Variables expected in the workspace:
%   sbml_path  : path to the SBML file to load
%   reps       : number of timed iterations
%   stop_time  : integration end time (days)
%   out_csv    : per-iteration timings CSV (one column, `seconds`)
%   load_csv   : path to write a single-row CSV with the model-load time
%
% MATLAB process startup is *not* timed here; that is measured separately
% by the orchestrating Python script via `matlab -batch` wall-clock.

required = {'sbml_path', 'reps', 'stop_time', 'out_csv', 'load_csv'};
for i = 1:numel(required)
    if ~exist(required{i}, 'var')
        error('bench_matlab:missingVar', '%s not set', required{i});
    end
end

% Optional: build a 6-bolus dose schedule that mirrors the C++ side's
% scenario.yaml. Same start/interval/count, same target species, same
% per-dose amount in mole.
if ~exist('apply_doses', 'var') || isempty(apply_doses)
    apply_doses = false;
end

t_load_start = tic;
model = sbmlimport(sbml_path);
cfg = getconfigset(model);
set(cfg, 'SolverType', 'sundials');
cfg.SolverOptions.RelativeTolerance = 1e-6;
cfg.SolverOptions.AbsoluteTolerance = 1e-9;
cfg.StopTime = stop_time;
cfg.SolverOptions.OutputTimes = (0:0.1:stop_time)';
t_load = toc(t_load_start);

dose_arg = {};
if apply_doses
    d = sbiodose('benchmark_dose', 'repeat');
    d.TargetName    = 'c1.Drug';
    d.StartTime     = 30;
    d.TimeUnits     = 'day';
    d.Amount        = 50;
    d.Interval      = 30;
    d.RepeatCount   = 5;       % 5 repeats after first dose => 6 total
    dose_arg = {d};
end

% Warm-up: first sbiosimulate triggers JIT/accelerator setup that we don't
% want to attribute to the integrator.
sbiosimulate(model, dose_arg{:});

times = zeros(reps, 1);
last_simdata = [];
for i = 1:reps
    t = tic;
    last_simdata = sbiosimulate(model, dose_arg{:});
    times(i) = toc(t);
end

% Optional: write the last trajectory to disk so the parity harness can
% compare it against the C++ side. Column names use the
% "compartment_species" form that compare_trajectories.py normalizes the
% C++ side ("compartment.species") to.
if exist('traj_csv', 'var') && ~isempty(traj_csv)
    % Bare species names repeat across compartments (every compartment
    % holds a species literally called "Drug"), so containers.Map collapses
    % them. Walk model.Species and last_simdata.DataNames in lockstep,
    % assuming simdata preserves model order — true for sbiosimulate.
    raw_names = last_simdata.DataNames;
    sp_list = model.Species;
    if length(raw_names) ~= length(sp_list)
        error('bench_matlab:simdataLen', ...
              'simdata column count (%d) != model species count (%d)', ...
              length(raw_names), length(sp_list));
    end
    col_names = cell(1, length(raw_names));
    for k = 1:length(raw_names)
        s = sp_list(k);
        if strcmp(raw_names{k}, s.Name)
            col_names{k} = [s.Parent.Name '_' s.Name];
        else
            col_names{k} = raw_names{k};
        end
    end
    T = array2table([last_simdata.Time, last_simdata.Data], ...
                    'VariableNames', [{'Time'}, col_names]);
    writetable(T, traj_csv);
end

writematrix(times, out_csv);
writematrix(t_load, load_csv);

fprintf('MATLAB load: %.4fs; integration median: %.4fs (n=%d)\n', ...
        t_load, median(times), reps);
