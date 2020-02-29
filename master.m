%% TWO-ASSET HOUSEHOLD MODEL
% This is the main script for this code repository
% HA model with a liquid asset and an illiquid asset
%
% Prior to running this script:
%
% (1) Set options in the section below. Directory only needs to be set for
% either the server or the local directory, depending on whether
% runopts.Server = 0 or 1
%
% (2) Modify the parameters script 'code/+setup/+params/get_params.m' and make sure that 
% runopts.mode is equal to 'get_params'. Note that all parameter defaults
% are set in the file HACTLib/+model_objects/ParamsDefaults.m, and get_params.m overrides
% these defaults. Any parameters not in get_params.m are set to their
% defaults. See the attributes of Params.m for a list of all
% parameters.
%
% (3) Set runopts.param_index equal to the index of the parameterization
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

runopts.calibrate = true;
runopts.Server = 1; % sets fast=0, param_index=slurm env var
runopts.fast = 0; % use small grid for debugging
runopts.mode = 'SDU_tests_new'; % 'get_params', 'grid_tests', 'chi0_tests', 'chi1_chi2_tests', 'table_tests', 'SDU_tests'
runopts.ComputeMPCS = true;
runopts.ComputeMPCS_illiquid = true;
runopts.SimulateMPCS = false; % also estimate MPCs by simulation
runopts.ComputeMPCS_news = false; % MPCs out of news, requires ComputeMPCS = 1
runopts.SimulateMPCS_news = false; % NOT CODED?

% whether or not to account for b = bmin, a > 0 case where household
% withdraws only enough to consume
runopts.DealWithSpecialCase = 0;

% Select which parameterization to run from parameters file
% (ignored when runops.Server = 1)
runopts.param_index = 1;

runopts.serverdir = '/home/livingstonb/GitHub/Continuous_Time_HA/';
runopts.localdir = '/home/brian/Documents/GitHub/Continuous_Time_HA/';


%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE
% -------------------------------------------------------------------------

if runopts.Server == 0
	runopts.direc = runopts.localdir;
else
	runopts.direc = runopts.serverdir;
	runopts.param_index = str2num(getenv('SLURM_ARRAY_TASK_ID'));
    
	runopts.fast = 0;
end

% check that specified directories exist
if ~exist(runopts.direc,'dir')
    error([runopts.direc ' does not exist'])
end
    
% for saving
runopts.suffix = num2str(runopts.param_index);

% directory to save output, and temp directory
runopts.savedir = fullfile(runopts.direc, 'output');

% temp directory
runopts.temp = fullfile(runopts.direc, 'temp');

addpath(fullfile(runopts.direc, 'code'));
addpath(fullfile(runopts.direc, 'factorization_lib'));

mkdir(runopts.temp);
mkdir(runopts.savedir);
addpath(runopts.temp);
addpath(runopts.savedir);

cd(runopts.direc)

%% --------------------------------------------------------------------
% GET PARAMETERS
% ---------------------------------------------------------------------
p = setup.params.(runopts.mode)(runopts);
p.print();

%% ------------------------------------------------------------------------
% CALIBRATING WITH SOLVER
% -------------------------------------------------------------------------
if ~isempty(p.calibrator)
	lbounds = p.calibrator.lbounds;
	ubounds = p.calibrator.ubounds;
    x0 = p.calibrator.x0;
    [calibrated_params, resnorm] = lsqnonlin(p.calibrator.solver_handle,...
    	x0, lbounds, ubounds);
    
    if resnorm > 1e-5
        error('Could not match targets')
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
stats = main(runopts, p);
toc

experiment.p = p;
experiment.stats = stats;
table_gen = HACTLib.model_objects.TableGenerator();
results_table = table_gen.create(experiment)
