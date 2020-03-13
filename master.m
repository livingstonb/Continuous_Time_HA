%% TWO-ASSET HOUSEHOLD MODEL
% This is the main script for this code repository
% HA model with a liquid asset and an illiquid asset
%
% Prior to running this script:
%
% (1) Set options in the section below. Directory only needs to be set for
% either the server or the local directory, depending on whether
% run_opts.Server = 0 or 1
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
% RUNNING ON THE SERVER: To run in batch on the server, use 
% code/batch/server.sbatch as a template. That script sends an array to SLURM 
% that runs all of the requested parameterizations in get_params.m. Output files
% are stored in the Output directory

clear all
warning('off', 'MATLAB:nearlySingularMatrix')

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------

param_opts.calibrate = true;
param_opts.fast = false; % use small grid for debugging
param_opts.ComputeMPCS = true;
param_opts.ComputeMPCS_illiquid = true; 
param_opts.SimulateMPCS = false; % also estimate MPCs by simulation
param_opts.ComputeMPCS_news = false;
param_opts.SimulateMPCS_news = false;
param_opts.DealWithSpecialCase = false;
param_opts.param_index = 1;

run_opts.Server = true;
run_opts.param_script = 'params_adj_cost_tests';
run_opts.serverdir = '/home/livingstonb/GitHub/Continuous_Time_HA/';
run_opts.localdir = '/home/brian/Documents/GitHub/Continuous_Time_HA/';

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE
% -------------------------------------------------------------------------

if ~run_opts.Server
	param_opts.direc = run_opts.localdir;
else
	param_opts.direc = run_opts.serverdir;
	param_opts.param_index = str2num(getenv('SLURM_ARRAY_TASK_ID'));
	param_opts.fast = false;
end

% temp directory
param_opts.temp_dir = fullfile(param_opts.direc, 'temp');

% output directory
param_opts.out_dir = fullfile(param_opts.direc, 'output');

fname = sprintf('output_%d.mat', param_opts.param_index);
param_opts.save_path = fullfile(param_opts.out_dir, fname);

addpath(fullfile(param_opts.direc, 'code'));
addpath(fullfile(param_opts.direc, 'factorization_lib'));

% check that specified directories exist
if ~exist(param_opts.direc, 'dir')
    error(strcat(param_opts.direc, ' does not exist'))
end

if ~exist(param_opts.temp_dir, 'dir')
    mkdir(param_opts.temp_dir);
end

if ~exist(param_opts.out_dir, 'dir')
    mkdir(param_opts.out_dir);
end

cd(param_opts.direc)
addpath(param_opts.temp_dir);
addpath(param_opts.out_dir);

%% --------------------------------------------------------------------
% GET PARAMETERS
% ---------------------------------------------------------------------
p = setup.params.(run_opts.param_script)(param_opts, param_opts.param_index);

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
	    options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1);
	    [calibrated_params, resnorm] = ...
	    	lsqnonlin(p.calibrator.solver_handle, x0, lbounds, ubounds, options);

	   	i_x0 = i_x0 + 1;
	end

    if (resnorm >= 1e-4)
        warning('Could not match targets')
    end
end

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION FILE
% -------------------------------------------------------------------------

%% r_b, r_a calibration
% rho calibrated to RA = 1
% if p.chi1 == 0.15
% 	[r_b_0, r_a_0] = setup.initial_returns_orig_adjcosts_ra1_calibration(p);
% else
% 	[r_b_0, r_a_0] = setup.initial_returns_new_adjcosts_ra1_calibration(p);
% end
% p.set("r_b", r_b_0);
% p.set("r_a", r_a_0);


%% returns risk, rho calibrated to RA = 5
% if p.chi1 == 0.15
% 	[r_b_0, r_a_0] = setup.initial_returns_orig_adjcosts_ra5_calibration(p);
% else
% 	[r_b_0, r_a_0] = setup.initial_returns_new_adjcosts_ra5_calibration(p);
% end


% final run
tic
save_results = true;
stats = main(p, save_results);
toc

% table_gen = HACTLib.tables.TableGenDetailed(p, stats);
% results_table = table_gen.create(p, stats)
table_gen = HACTLib.tables.TableFancy(p, {stats});
results_table = table_gen.create(p, {stats})

xlx_path = sprintf('run%d_table.xlsx', p.param_index);
xlx_path = fullfile(p.out_dir, xlx_path);
writetable(results_table, xlx_path, 'WriteRowNames', true)
