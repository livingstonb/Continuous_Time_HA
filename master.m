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
warning('off','MATLAB:nearlySingularMatrix')

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------

runopts.calibrate = true;
runopts.Server = 1; % sets IterateRho=1,fast=0,param_index=slurm env var
runopts.fast = 0; % use small grid for debugging
runopts.mode = 'params_adj_cost_tests'; % 'get_params', 'grid_tests', 'chi0_tests', 'chi1_chi2_tests', 'table_tests', 'SDU_tests'
runopts.ComputeMPCS = true;
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
if strcmp(runopts.mode,'grid_test')
	p = setup.params.grid_test_params(runopts);
elseif strcmp(runopts.mode,'chi0_tests')
    p = setup.params.chi0_tests(runopts);
elseif strcmp(runopts.mode,'chi1_chi2_tests')
    p = setup.params.chi1_chi2_tests(runopts);
elseif strcmp(runopts.mode,'table_tests')
    p = setup.params.table_tests(runopts);
elseif strcmp(runopts.mode,'table_tests_bequests')
    p = setup.params.table_tests_bequests(runopts);
elseif strcmp(runopts.mode,'get_params')
	p = setup.params.get_params(runopts);
elseif strcmp(runopts.mode, 'SDU_tests')
    p = setup.params.SDU_tests(runopts);
elseif strcmp(runopts.mode, 'endog_labor_tests')
	p = setup.params.endog_labor_tests(runopts);
elseif strcmp(runopts.mode, 'params_one_asset')
	p = setup.params.params_one_asset(runopts);
elseif strcmp(runopts.mode, 'params_adj_cost_tests')
	p = setup.params.params_adj_cost_tests(runopts);
end
p.print();

%% ------------------------------------------------------------------------
% CALIBRATING WITH FSOLVE
% -------------------------------------------------------------------------
if ~isempty(p.calibrator)
    calibrated_params = fsolve(p.calibrator, p.x0_calibration);
end

% %% ------------------------------------------------------------------------
% % CALIBRATING WITH FSOLVE
% % -------------------------------------------------------------------------
% import HACTLib.model_objects.AltCalibrator
% n_calibrations = 0;

% % Vary (rho, r_a) to match median(a+b) = 1.6, median(b) = 0.1
% param_name = {'rho', 'r_a'};
% stat_name = {'median_totw', 'median_liqw'};
% stat_target = [1.6, 0.1];
% inits = [p.rho, p.r_a];
% n_calibrations = n_calibrations + 1;

% if p.OneAsset
% 	param_name = param_name(1);
% 	stat_name = stat_name(1);
% 	stat_target = stat_target(1);
% 	inits = inits(1);
% end

% rho_bounds = [0.002, 0.4];
% ra_bounds = [p.r_b+1e-4, 0.3];

% if n_calibrations == 1
%     calibrator = AltCalibrator(p, runopts, param_name,...
%         stat_name, stat_target);

%     calibrator.set_param_bounds(1, rho_bounds);
%     if ~p.OneAsset
%     	calibrator.set_param_bounds(2, ra_bounds);
%     end

%     for ii = 1:numel(inits)
%     	z_inits(ii) = calibrator.convert_to_solver_input(inits(ii), ii);
%     end

%     calibrated_params = fsolve(@(x) calibrator.fn_handle(x, p), z_inits(:));
% elseif n_calibrations > 1
%     error("Ensure that a max of one calibration is selected")
% end

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
% 
% % start iteration
% calibrator = HACTLib.computation.Calibrator(runopts, p, "r_b, r_a");
% x0 = calibrator.create_initial_condition([r_b_0, r_a_0]);
% 
% opts = optimoptions('fsolve', 'MaxFunctionEvaluations', 400, 'MaxIterations', 600);
% fsolve(calibrator.objective, x0, opts);
% p.set("NoRisk", 1);
% p.set("ComputeMPCS", 1);
% p.set("SaveResults", 1);


%% r_a, rho calibration
% calibrator = HACTLib.computation.Calibrator(runopts, p, "r_a, rho");
% x0 = calibrator.create_initial_condition([0.066513, 0.142622]);
% fsolve(calibrator.objective, x0);
% p.set("NoRisk", 1);
% p.set("ComputeMPCS", 1);
% p.set("SaveResults", 1);

%% rho calibration
% calibrator = HACTLib.computation.Calibrator(runopts, p, "rho");
% x0 = calibrator.create_initial_condition(0.035);
% fsolve(calibrator.objective, x0);
% p.set("NoRisk", 1);
% p.set("ComputeMPCS", 1);
% p.set("SaveResults", 1);

% final run
tic
stats = main(runopts, p);
toc

experiment.p = p;
experiment.stats = stats;
table_gen = HACTLib.model_objects.TableGenerator();
results_table = table_gen.create(experiment)
