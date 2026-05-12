% export_matlab_trajectories.m - Run SimBiology model and export trajectories
%
% Usage:
%   matlab -batch "output_csv='/path/to/out.csv'; matlab_model_dir='/path/to/repo'; matlab_model_script='model_script_name'; run('/path/to/export_matlab_trajectories.m')"
%
% `matlab_model_dir` is the consumer repo's root (the startup.m path). The
% `pdac_build_dir` name is accepted for backward compatibility.
%
% Requires: output_csv variable set before calling

if ~exist('output_csv', 'var')
    error('Set output_csv before running this script');
end
% matlab_model_dir is the preferred name; fall back to the legacy
% pdac_build_dir so existing callers keep working during migration.
if ~exist('matlab_model_dir', 'var') || isempty(matlab_model_dir)
    if exist('pdac_build_dir', 'var') && ~isempty(pdac_build_dir)
        matlab_model_dir = pdac_build_dir;
    else
        matlab_model_dir = pwd;
    end
end

cd(matlab_model_dir);
run('startup.m');
% Add this test dir to path so yaml_read.m is callable.
addpath(fileparts(mfilename('fullpath')));

% Use the live .m model script. sbmlimport(sbml_path) would be the
% source-of-truth loader (same SBML the C++ codegen reads), but SimBiology's
% sbmlexport strips the max(x, 0) guards from sqrt/pow kinetic laws in
% certain kinetic expressions — reimported models can produce complex-valued
% RHS and fail to integrate. Live .m it is. The caller sets
% `matlab_model_script` to the bare script name (no `.m`); it is `run`-ed so
% its workspace (notably `model`) lands in ours.
if ~exist('matlab_model_script', 'var') || isempty(matlab_model_script)
    error(['Set matlab_model_script before running this script ' ...
           '(e.g. matlab_model_script = ''immune_oncology_model_PDAC'').']);
end
run(matlab_model_script);

% Tighten solver tolerances to roughly match the C++ CVODE settings
% (reltol=1e-6, abstol=1e-12). SimBiology defaults (1e-3 / 1e-6) are
% too loose to compare against C++ at day-100+ timescales.
cfg = getconfigset(model);
% sundials is what pdac-build uses in production; the default (ode15s) chokes
% on the stiff bolus transient in some scenarios.
set(cfg, 'SolverType', 'sundials');
cfg.SolverOptions.RelativeTolerance = 1e-6;
cfg.SolverOptions.AbsoluteTolerance = 1e-9;
% Output grid. Defaults match the baseline 365-day run; set `stop_time`
% before calling to use a shorter scenario (e.g. event-fire tests).
if ~exist('stop_time', 'var') || isempty(stop_time)
    stop_time = 365;
end
cfg.StopTime = stop_time;
% If output_times_csv is provided, pin MATLAB's sample grid to the C++
% grid so the parity comparison is row-aligned (no interpolation across
% bolus discontinuities). Otherwise use the default 0.1 d grid.
if exist('output_times_csv', 'var') && ~isempty(output_times_csv)
    cfg.SolverOptions.OutputTimes = readmatrix(output_times_csv);
else
    cfg.SolverOptions.OutputTimes = (0:0.1:stop_time)';
end

% Optional: override model values with those from a param_all XML so MATLAB
% simulates with the exact same ICs / parameters as the C++ run.
if exist('param_xml', 'var') && ~isempty(param_xml)
    fprintf('Applying param overrides from %s\n', param_xml);
    doc = xmlread(param_xml);
    iv = doc.getElementsByTagName('init_value').item(0);
    if isempty(iv)
        error('No <init_value> element found in %s', param_xml);
    end

    % Index model objects by the names we expect in the XML
    comp_map = containers.Map();
    for i = 1:length(model.Compartments)
        comp_map(model.Compartments(i).Name) = model.Compartments(i);
    end
    sp_map = containers.Map();
    for i = 1:length(model.Species)
        s = model.Species(i);
        sp_map([s.Parent.Name '_' s.Name]) = s;
    end
    par_map = containers.Map();
    for i = 1:length(model.Parameters)
        par_map(model.Parameters(i).Name) = model.Parameters(i);
    end

    n_set = struct('comp', 0, 'sp', 0, 'par', 0, 'miss', 0);
    for section_name = {'Compartment', 'Species', 'Parameter'}
        sec = iv.getElementsByTagName(section_name{1}).item(0);
        if isempty(sec); continue; end
        children = sec.getChildNodes();
        for k = 0:children.getLength()-1
            node = children.item(k);
            if node.getNodeType() ~= node.ELEMENT_NODE; continue; end
            name = char(node.getNodeName());
            val = str2double(char(node.getTextContent()));
            if isnan(val); continue; end
            switch section_name{1}
                case 'Compartment'
                    if comp_map.isKey(name)
                        c = comp_map(name); c.Capacity = val;
                        n_set.comp = n_set.comp + 1;
                    else; n_set.miss = n_set.miss + 1; end
                case 'Species'
                    if sp_map.isKey(name)
                        s = sp_map(name); s.InitialAmount = val;
                        n_set.sp = n_set.sp + 1;
                    else; n_set.miss = n_set.miss + 1; end
                case 'Parameter'
                    if par_map.isKey(name)
                        p = par_map(name); p.Value = val;
                        n_set.par = n_set.par + 1;
                    else; n_set.miss = n_set.miss + 1; end
            end
        end
    end
    fprintf('  Overrode: %d compartments, %d species, %d parameters (%d names not in model)\n', ...
        n_set.comp, n_set.sp, n_set.par, n_set.miss);
end

% Optional: dosing scenario. If scenario_yaml is set, build a SimBiology
% dose_schedule via schedule_dosing.m using the scenario's drug list +
% overrides, and pass it to sbiosimulate. drugs/dose/schedule/patient
% values are in the same format schedule_dosing.m expects.
dose_schedule = [];
if exist('scenario_yaml', 'var') && ~isempty(scenario_yaml)
    fprintf('Loading dosing scenario from %s\n', scenario_yaml);
    scenario = yaml_read(scenario_yaml);
    % Normalize drug list before the empty check: yaml_read can parse a
    % YAML `drugs: []` as `{[]}` (a 1×1 cell containing []), which the
    % cheap `~isempty(...)` guard misses, sending an empty drug name into
    % schedule_dosing and crashing drugless scenarios.
    drug_list = {};
    if isfield(scenario, 'dosing') && isfield(scenario.dosing, 'drugs')
        drug_list = scenario.dosing.drugs;
        if ~iscell(drug_list); drug_list = {drug_list}; end
        drug_list = drug_list(~cellfun(@(d) ...
            isempty(d) || (ischar(d) && isempty(strtrim(d))), drug_list));
    end
    if isempty(drug_list)
        fprintf('  no drugs scheduled\n');
    else
        dosing_args = {};
        patient_weight = 70;
        patient_bsa = 1.9;
        if isfield(scenario.dosing, 'patientWeight')
            patient_weight = scenario.dosing.patientWeight;
        end
        if isfield(scenario.dosing, 'patientBSA')
            patient_bsa = scenario.dosing.patientBSA;
        end

        % Pass through each <drug>_dose and <drug>_schedule pair.
        % Schedule arrays in YAML are [start, interval, num_doses];
        % schedule_dosing.m wants [start, interval, RepeatCount] where
        % RepeatCount = num_doses - 1.
        dose_fields = fieldnames(scenario.dosing);
        for ii = 1:numel(dose_fields)
            f = dose_fields{ii};
            if any(strcmp(f, {'drugs', 'patientWeight', 'patientBSA'})); continue; end
            val = scenario.dosing.(f);
            if endsWith(f, '_schedule')
                % Unpack cell → numeric triple.
                if iscell(val)
                    sched = cellfun(@(x) x, val);
                else
                    sched = val;
                end
                if numel(sched) ~= 3
                    error('%s must be [start, interval, num_doses]', f);
                end
                sched(3) = max(0, sched(3) - 1);  % num_doses → RepeatCount
                val = sched;
            end
            dosing_args{end+1} = f;
            dosing_args{end+1} = val;
        end
        dosing_args{end+1} = 'patientWeight';
        dosing_args{end+1} = patient_weight;
        dosing_args{end+1} = 'patientBSA';
        dosing_args{end+1} = patient_bsa;
        dose_schedule = schedule_dosing(drug_list, dosing_args{:});
        fprintf('  %d dose object(s) scheduled\n', numel(dose_schedule));
    end
end

% Optional: run a consumer-provided natural-history / IC-evolution
% function (YAML-driven; reads the same healthy-state file the C++ side
% consumes). Replaces the model ICs with the diagnosis-time state;
% the scenario then runs on top of that. Mirrors the C++ dumper's
% --evolve-to-diagnosis flow so both sides solve the same problem.
%
% Consumers set `matlab_evolve_function` to a function handle with the
% signature [model, ok, info] = fn(model, varargin). Backward-compat:
% if `evolve_to_diagnosis_enabled` is set but no handle is provided,
% call the legacy `evolve_to_diagnosis` name on MATLAB's path.
if exist('matlab_evolve_function', 'var') && ~isempty(matlab_evolve_function)
    evolve_fn = matlab_evolve_function;
    evolve_enabled = true;
elseif exist('evolve_to_diagnosis_enabled', 'var') && evolve_to_diagnosis_enabled
    evolve_fn = @evolve_to_diagnosis;
    evolve_enabled = true;
else
    evolve_enabled = false;
end
if evolve_enabled
    fprintf('Running consumer evolve function...\n');
    [model, evolve_ok, ~] = evolve_fn(model, 'Debug', true);
    if ~evolve_ok
        error('evolve function rejected this parameter set');
    end
    % Reinstate the scenario's output grid in case the evolve mutated it.
    cfg.StopTime = stop_time;
    if exist('output_times_csv', 'var') && ~isempty(output_times_csv)
        cfg.SolverOptions.OutputTimes = readmatrix(output_times_csv);
    else
        cfg.SolverOptions.OutputTimes = (0:0.1:stop_time)';
    end
end

if isempty(dose_schedule)
    simdata = sbiosimulate(model);
else
    simdata = sbiosimulate(model, dose_schedule);
end

t = simdata.Time;
data = simdata.Data;
raw_names = simdata.DataNames;

% Map bare species name -> qualified "Compartment_Species". A species name may
% repeat across compartments, so keep a list per bare name and index into it.
species_list = model.Species;
bare_to_comp = containers.Map();

for i = 1:length(species_list)
    s = species_list(i);
    comp = s.Parent.Name;
    qname = [comp '_' s.Name];
    if bare_to_comp.isKey(s.Name)
        existing = bare_to_comp(s.Name);
        if ~iscell(existing)
            existing = {existing};
        end
        existing{end+1} = qname;
        bare_to_comp(s.Name) = existing;
    else
        bare_to_comp(s.Name) = qname;
    end
end

n_cols = length(raw_names);
col_names = cell(1, n_cols);
seen_count = containers.Map();

for i = 1:n_cols
    n = raw_names{i};
    if bare_to_comp.isKey(n)
        mapping = bare_to_comp(n);
        if iscell(mapping)
            if seen_count.isKey(n)
                seen_count(n) = seen_count(n) + 1;
            else
                seen_count(n) = 1;
            end
            idx = seen_count(n);
            if idx <= length(mapping)
                col_names{i} = mapping{idx};
            else
                col_names{i} = [n '_' num2str(idx)];
            end
        else
            if seen_count.isKey(n)
                col_names{i} = n;
            else
                seen_count(n) = 1;
                col_names{i} = mapping;
            end
        end
    else
        col_names{i} = n;
    end
end

final_seen = containers.Map();
for i = 1:n_cols
    n = col_names{i};
    if final_seen.isKey(n)
        final_seen(n) = final_seen(n) + 1;
        col_names{i} = [n '_dup' num2str(final_seen(n))];
    else
        final_seen(n) = 1;
    end
end

T = array2table([t, data], 'VariableNames', [{'Time'}, col_names]);
writetable(T, output_csv);
fprintf('Wrote %d time points, %d columns to %s\n', size(data,1), n_cols, output_csv);
