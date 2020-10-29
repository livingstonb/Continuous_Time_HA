%% TWO-ASSET HOUSEHOLD MODEL
% This is the main script for this code repository
% HA model with a liquid asset and an illiquid asset
%
% Prior to running this script:
%
% (1) Set options in the section below.
%
% (2) Modify the parameters script 'code/+setup/+params/get_params.m' and make sure that 
% run_opts.mode is equal to 'get_params'. Note that all parameter defaults
% are set in the file HACTLib/+model_objects/ParamsDefaults.m, and get_params.m overrides
% these defaults. Any parameters not in get_params.m are set to their
% defaults. See the attributes of Params.m for a list of all
% parameters.
%
% (3) Set run_opts.param_index equal to the index of the parameterization
% you would like to run in the parameters file (1,2,...).
%
% (4) cd into the Continuous_Time_HA directory.
%
% RUNNING ON THE SERVER: To run in batch on the server, use 
% code/batch/server.sbatch as a template. That script sends an array to SLURM 
% that runs all of the requested parameterizations in get_params.m. Output files
% are stored in the Output directory

clear all
warning('off', 'MATLAB:nearlySingularMatrix')

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------

param_opts.calibrate = false;
param_opts.fast = true; % use small grid for debugging
param_opts.ComputeMPCS = true;
param_opts.ComputeMPCS_illiquid = false; 
param_opts.SimulateMPCS = false; % also estimate MPCs by simulation
param_opts.ComputeMPCS_news = false;
param_opts.SimulateMPCS_news = false;
param_opts.DealWithSpecialCase = false; % need to recode this
param_opts.param_index = 1;
param_opts.makePlots = false; % not coded yet

run_opts.check_nparams = false;
run_opts.Server = false;
run_opts.param_script = 'params_adj_cost_tests';

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE
% -------------------------------------------------------------------------
[~, currdir] = fileparts(pwd());
if ~strcmp(currdir, 'Continuous_Time_HA')
    msg = 'The user must cd into the Continuous_Time_HA directory';
    bad_dir = MException('Continuous_Time_HA:master', msg);
    throw(bad_dir);
end

if run_opts.Server
	param_opts.param_index = str2num(getenv('SLURM_ARRAY_TASK_ID'));
	param_opts.fast = false;
    run_opts.check_nparams = false;
end

addpath('code');
addpath('factorization_lib');

warning('off', 'MATLAB:MKDIR:DirectoryExists')
mkdir('temp');
mkdir('output');

%% --------------------------------------------------------------------
% GET PARAMETERS
% ---------------------------------------------------------------------
[p, nparams] = setup.params.(run_opts.param_script)(param_opts);
if run_opts.check_nparams
    fprintf('Parameters script contains %d specifications\n', nparams)
    return
end

% Create Params object
p = HACTLib.model_objects.Params(p);
p.print();

%% ------------------------------------------------------------------------
% CALIBRATING WITH SOLVER
% -------------------------------------------------------------------------
if ~isempty(p.calibrator)
	lbounds = p.calibrator.lbounds;
	ubounds = p.calibrator.ubounds;

	resnorm = 100;
	n_x0 = numel(p.calibrator.x0);
	i_x0 = 1;
	while (resnorm >= 1e-4) && (i_x0 <= n_x0)
	    x0 = p.calibrator.x0{i_x0};
	    options = optimoptions(@lsqnonlin,...
	    	'MaxFunctionEvaluations', p.calibration_maxiters,...
	    	'FunctionTolerance', p.calibration_crit,...
	    	'OptimalityTolerance', p.calibration_crit);
	    [calibrated_params, resnorm] = ...
	    	lsqnonlin(p.calibrator.solver_handle, x0, lbounds, ubounds, options);

	   	i_x0 = i_x0 + 1;
	end

    if (resnorm >= 1e-4)
        error('Could not match targets')
    end
end

save_results = true;
stats = main(p, save_results);

table_gen = HACTLib.tables.StatsTable(p, {stats});
results_table = table_gen.create(p, {stats})

xlx_path = sprintf('run%d_table.xlsx', p.param_index);
xlx_path = fullfile('output', xlx_path);
writetable(results_table, xlx_path, 'WriteRowNames', true)
