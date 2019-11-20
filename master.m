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
% (2) Modify the parameters script 'get_params.m' and make sure that 
% runopts.mode is equal to 'get_params'. Note that all parameter defaults
% are set in the class file code/+setup/Params.m, and get_params.m overrides
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

runopts.Server = 0; % sets IterateRho=1,fast=0,param_index=slurm env var
runopts.fast = 1; % use small grid for debugging
runopts.mode = 'get_params'; % 'get_params', 'grid_tests', 'chi0_tests', 'chi1_chi2_tests', 'table_tests', 'SDU_tests'
runopts.ComputeMPCS = 0;
runopts.SimulateMPCS = 0; % also estimate MPCs by simulation
runopts.ComputeMPCS_news = 0; % MPCs out of news, requires ComputeMPCS = 1
runopts.SimulateMPCS_news = 0; % NOT CODED

% whether or not to account for b = bmin, a > 0 case where household
% withdraws only enough to consume
runopts.DealWithSpecialCase = 0;  

% Select which parameterization to run from parameters file
% (ignored when runops.Server = 1)
runopts.param_index = 2;

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
runopts.savedir = [runopts.direc 'output/'];

% temp directory
runopts.temp = [runopts.direc 'temp/'];

addpath([runopts.direc,'code']);
addpath([runopts.direc,'code/factorization_lib']);

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
end
p.print();

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION FILE
% -------------------------------------------------------------------------


% % to calibrate to (rb, ra)
% calibrator = solver.Calibrators.rb_ra_calibrator(runopts, p);
% if p.invies == 1
%     switch p.riskaver
%         case 1
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [0.005, 0.022]);
%         case 2
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [-0.04, 0.0001]);
%         case 5
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [-0.06, -0.001]);
%         case 10
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [-0.1, -0.002]);
%         case 20
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [-0.12, -0.005]);
%     end
% else
%     switch p.riskaver
%         case 1
%             p.set("riskaver", 1.01);
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [0.004, 0.012]);
%         case 2
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [-0.035, 0.0015]);
%         case 5
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [-0.085, -0.0003]);
%         case 10
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [-0.125, -0.001]);
%         case 20
%             x0 = solver.Calibrators.rb_ra_get_initial(p, [-0.15, -0.01]);
%     end
% end
% fsolve(calibrator, x0);

% calibrator = solver.Calibrators.ra_rho_calibrator(runopts, p);
% x0 = solver.Calibrators.ra_rho_get_initial(p, [0.033944, 0.022869]);
% fsolve(calibrator, x0);

calibrator = solver.Calibrators.rho_calibrator(runopts, p);
fsolve(calibrator, [0.8]);

% final run
p.set("NoRisk", 1);
p.set("ComputeMPCS", 1);
% p.set("ComputeMPCS_news", 1);
p.set("SaveResults", 1);
tic
stats = main(runopts, p);
toc
