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
runopts.IterateRho = 1; % if set to zero, the parameter 'rho' is used
runopts.fast = 0; % use small grid for debugging
runopts.mode = 'SDU_tests'; % 'get_params', 'grid_tests', 'chi0_tests', 'chi1_chi2_tests', 'table_tests', 'SDU_tests'
runopts.ComputeMPCS = 0;
runopts.SimulateMPCS = 0; % also estimate MPCs by simulation
runopts.ComputeMPCS_news = 0; % MPCs out of news, requires ComputeMPCS = 1
runopts.SimulateMPCS_news = 0; % NOT CODED

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
runopts.savedir = [runopts.direc 'output/two_asset/'];

% temp directory
runopts.temp = [runopts.direc 'temp/two_asset/'];

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
	p = setup.two_asset.params.grid_test_params(runopts);
elseif strcmp(runopts.mode,'chi0_tests')
    p = setup.two_asset.params.chi0_tests(runopts);
elseif strcmp(runopts.mode,'chi1_chi2_tests')
    p = setup.two_asset.params.chi1_chi2_tests(runopts);
elseif strcmp(runopts.mode,'table_tests')
    p = setup.two_asset.params.table_tests(runopts);
elseif strcmp(runopts.mode,'table_tests_bequests')
    p = setup.two_asset.params.table_tests_bequests(runopts);
elseif strcmp(runopts.mode,'get_params')
	p = setup.two_asset.params.get_params(runopts);
elseif strcmp(runopts.mode, 'SDU_tests')
    p = setup.two_asset.params.SDU_tests(runopts);
end
p.print();

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION FILE
% -------------------------------------------------------------------------


% % to calibrate to (rb, ra) (turn off rho iteration)
% calibrator = @(r) solver.two_asset.risk_premium_calibrator(r, runopts, p);
% % returns = fsolve(calibrator, log([0.02/4+0.05, 0.04/4]));
% returns = fsolve(calibrator, [0.4, 0.5]);
% 
% new_rb = 0.035*(returns(1))/(1+abs(returns(1)));
% new_ra = new_rb + 0.04 * abs(returns(2)) / (1 + abs(returns(2)));
% 
% p.reset_returns(new_rb, new_ra);
% fprintf("FINAL LIQUID RETURN = %f\n", p.r_b)
% fprintf("FINAL ILLIQUID RETURN = %f\n", p.r_a)
% 


% to calibrate to (rho, ra) (turn on rho iteration)
calibrator = @(ra) solver.two_asset.ra_calibrator(ra, runopts, p);
returns = fsolve(calibrator, log(0.06/4));
p.reset_returns(p.r_b, exp(returns));
fprintf("FINAL RHO = %f\n", p.rho);
fprintf("FINAL ILLIQUID RETURN = %f\n", p.r_a);

% % final run
% tic
% stats = main_two_asset(runopts, p);
% toc
