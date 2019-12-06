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
runopts.savedir = fullfile(runopts.direc, 'output');

% temp directory
runopts.temp = fullfile(runopts.direc, 'temp');

addpath(fullfile(runopts.direc, 'code'));
addpath(fullfile(runopts.direc, 'code', 'factorization_lib'));

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

%% r_b, r_a calibration for rho calibrated to RA = 1, original adj costs

% if p.chi1 == 0.15
%     
% if p.invies == 1
%     switch p.riskaver
%         case 1
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = 0.005;
%                     r_a_0 = 0.022866;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = 0.005006;
%                     r_a_0 = 0.022877;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = 0.005026;
%                     r_a_0 = 0.022909;
%                 case 0.05
%                     converged = true;
%                     r_b_0 = 0.005154;
%                     r_a_0 = 0.02313;
%                 case 0.1
%                     converged = true;
%                     r_b_0 = 0.005539;
%                     r_a_0 = 0.023891;
%                 case 0.15
%                     converged = true;
%                     r_b_0 = 0.006049;
%                     r_a_0 = 0.025091;
%             end
%         case 2
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.009868;
%                     r_a_0 = 0.01987;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.009856;
%                     r_a_0 = 0.019937;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = -0.00982;
%                     r_a_0 = 0.020134;
%                 case 0.05
%                     r_b_0 = -0.0098;
%                     r_a_0 = 0.022134;
%                 case 0.1
%                     converged = true;
%                     r_b_0 = -0.010709;
%                     r_a_0 = 0.024941;
%                 case 0.15
%                     r_b_0 = -0.00961;
%                     r_a_0 = 0.032941;
%             end
%         case 5
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.051442;
%                     r_a_0 = 0.009666;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.051958;
%                     r_a_0 = 0.009898;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = -0.05323;
%                     r_a_0 = 0.010561;
%                 case 0.05
%                     converged = true;
%                     r_b_0 = -0.057679;
%                     r_a_0 = 0.014174;
%                 case 0.1
%                     converged = true;
%                     r_b_0 = -0.065535;
%                     r_a_0 = 0.022612;
%                 case 0.15
%                     converged = true;
%                     r_b_0 = -0.07403;
%                     r_a_0 = 0.03228;
%             end
%         case 10
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.092879;
%                     r_a_0 = 0.00526;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.094133;
%                     r_a_0 = 0.005732;
%                 case 0.02
%                     r_b_0 = -0.098433;
%                     r_a_0 = 0.009332;
%                 case 0.05
%                     r_b_0 = -0.105433;
%                     r_a_0 = 0.015832;
%                 case 0.1
%                     converged = true;
%                     r_b_0 = -0.118208;
%                     r_a_0 = 0.02589;
%                 case 0.15
%                     r_b_0 = -0.131433;
%                     r_a_0 = 0.036032;
%             end
%         case 20
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.124384;
%                     r_a_0 = 0.004371;
%                 case 0.01
%                     r_b_0 = -0.128384;
%                     r_a_0 = 0.006371;
%                 case 0.02
%                     r_b_0 = -0.133384;
%                     r_a_0 = 0.012371;
%                 case 0.05
%                     r_b_0 = -0.144384;
%                     r_a_0 = 0.018371;
%                 case 0.1
%                     r_b_0 = -0.142;
%                     r_a_0 = 0.026471;
%                 case 0.15
%                     r_b_0 = -0.1689;
%                     r_a_0 = 0.038;
%             end
%     end
% else
%     switch p.riskaver
%         case 1.01
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = 0.005000;
%                     r_a_0 = 0.019368;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = 0.00503;
%                     r_a_0 = 0.019402;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = 0.005118;
%                     r_a_0 = 0.019501;
%                 case 0.05
%                     r_b_0 = 0.005154;
%                     r_a_0 = 0.0199;
%                 case 0.1
%                     r_b_0 = 0.005539;
%                     r_a_0 = 0.021891;
%                 case 0.15
%                     r_b_0 = 0.006049;
%                     r_a_0 = 0.023091;
%             end
%         case 2
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.009868;
%                     r_a_0 = 0.01987;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.009856;
%                     r_a_0 = 0.019937;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = -0.00982;
%                     r_a_0 = 0.020134;
%                 case 0.05
%                     r_b_0 = -0.0098;
%                     r_a_0 = 0.022134;
%                 case 0.1
%                     converged = true;
%                     r_b_0 = -0.010709;
%                     r_a_0 = 0.024941;
%                 case 0.15
%                     r_b_0 = -0.00961;
%                     r_a_0 = 0.032941;
%             end
%         case 5
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.051442;
%                     r_a_0 = 0.009666;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.051958;
%                     r_a_0 = 0.009898;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = -0.05323;
%                     r_a_0 = 0.010561;
%                 case 0.05
%                     converged = true;
%                     r_b_0 = -0.057679;
%                     r_a_0 = 0.014174;
%                 case 0.1
%                     converged = true;
%                     r_b_0 = -0.065535;
%                     r_a_0 = 0.022612;
%                 case 0.15
%                     converged = true;
%                     r_b_0 = -0.07403;
%                     r_a_0 = 0.03228;
%             end
%         case 10
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.08399;
%                     r_a_0 = 0.00631;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.084389;
%                     r_a_0 = 0.006856;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = -0.084691;
%                     r_a_0 = 0.006332;
%                 case 0.05
%                     r_b_0 = -0.085239;
%                     r_a_0 = 0.008287;
%                 case 0.1
%                     r_b_0 = -0.097537;
%                     r_a_0 = 0.031703;
%                 case 0.15
%                     r_b_0 = -0.105433;
%                     r_a_0 = 0.041032;
%             end
%         case 20
%             switch p.sigma_r
%                 case 0
%                     r_b_0 = -0.101910;
%                     r_a_0 = 0.007058;
%                 case 0.01
%                     r_b_0 = -0.103299;
%                     r_a_0 = 0.007687;
%                 case 0.02
%                     r_b_0 = -0.104012;
%                     r_a_0 = 0.010050;
%                 case 0.05
%                     r_b_0 = -0.107618;
%                     r_a_0 = 0.028151;
%                 case 0.1
%                     r_b_0 = -0.148;
%                     r_a_0 = 0.030471;
%                 case 0.15
%                     r_b_0 = -0.1689;
%                     r_a_0 = 0.042;
%             end
%     end
% end
% 
% calibrator = solver.Calibrator(runopts, p, "r_b, r_a");
% x0 = calibrator.create_initial_condition([r_b_0, r_a_0]);
% 
% opts = optimoptions('fsolve', 'MaxFunctionEvaluations', 400, 'MaxIterations', 600);
% fsolve(calibrator.objective, x0, opts);
% 
% end
% 
% %% r_b, r_a calibration for rho calibrated to RA = 1, NEW adj costs
% 
% if p.chi1 ~= 0.15
%     
% if p.invies == 1
%     switch p.riskaver
%         case 1
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = 0.005;
%                     r_a_0 = 0.022612;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = 0.005006;
%                     r_a_0 = 0.022627;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = 0.005026;
%                     r_a_0 = 0.020909;
%                 case 0.05
%                     converged = true;
%                     r_b_0 = 0.005054;
%                     r_a_0 = 0.022669;
%                 case 0.1
%                     converged = true;
%                     r_b_0 = 0.005339;
%                     r_a_0 = 0.023991;
%                 case 0.15
%                     converged = true;
%                     r_b_0 = 0.005408;
%                     r_a_0 = 0.025668;
%             end
%         case 2
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.006883;
%                     r_a_0 = 0.019505;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.006915;
%                     r_a_0 = 0.019600;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = -0.007010;
%                     r_a_0 = 0.019881;
%                 case 0.05
%                     converged = true;
%                     r_b_0 = -0.007469;
%                     r_a_0 = 0.021652;
%                 case 0.1
%                     converged = true;
%                     r_b_0 = -0.009050;
%                     r_a_0 = 0.026614;
%                 case 0.15
%                     converged = true;
%                     r_b_0 = -0.013651; % converged, error with no risk
%                     r_a_0 = 0.033504; % converged, error with no risk
%             end
%         case 5
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.040727;
%                     r_a_0 = 0.008936;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.041411;
%                     r_a_0 = 0.009220;
%                 case 0.02
%                     r_b_0 = -0.043032;
%                     r_a_0 = 0.011039;
%                 case 0.05
%                     r_b_0 = -0.047679;
%                     r_a_0 = 0.015174;
%                 case 0.1
%                     r_b_0 = -0.053535;
%                     r_a_0 = 0.021612;
%                 case 0.15
%                     r_b_0 = -0.063892;
%                     r_a_0 = 0.030143;
%             end
%         case 10
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.076174;
%                     r_a_0 = 0.003779;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.078166;
%                     r_a_0 = 0.004316;
%                 case 0.02
%                     r_b_0 = -0.082166;
%                     r_a_0 = 0.005116;
%                 case 0.05
%                     r_b_0 = -0.090433;
%                     r_a_0 = 0.006132;
%                 case 0.1
%                     r_b_0 = -0.099433;
%                     r_a_0 = 0.007032;
%                 case 0.15
%                     r_b_0 = -0.111433;
%                     r_a_0 = 0.008032;
%             end
%         case 20
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.104539;
%                     r_a_0 = 0.002535;
%                 case 0.01
%                     r_b_0 = -0.107539;
%                     r_a_0 = 0.003371;
%                 case 0.02
%                     r_b_0 = -0.111539;
%                     r_a_0 = 0.004071;
%                 case 0.05
%                     r_b_0 = -0.114384;
%                     r_a_0 = 0.005371;
%                 case 0.1
%                     r_b_0 = -0.1252;
%                     r_a_0 = 0.006971;
%                 case 0.15
%                     r_b_0 = -0.1289;
%                     r_a_0 = 0.012;
%             end
%     end
% else
%     switch p.riskaver
%         case 1.01
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = 0.005;
%                     r_a_0 = 0.019791;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = 0.005026;
%                     r_a_0 = 0.019836;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = 0.005104;
%                     r_a_0 = 0.019966;
%                 case 0.05
%                     converged = true;
%                     r_b_0 = 0.005612;
%                     r_a_0 = 0.020795;
%                 case 0.1
%                     r_b_0 = 0.0063308;
%                     r_a_0 = 0.023997;
%                 case 0.15
%                     r_b_0 = 0.0071308;
%                     r_a_0 = 0.030597;
%             end
%         case 2
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.005964;
%                     r_a_0 = 0.017394;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.005943;
%                     r_a_0 = 0.017522;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = -0.005874;
%                     r_a_0 = 0.017893;
%                 case 0.05
%                     converged = true;
%                     r_b_0 = -0.005546;
%                     r_a_0 = 0.020098;
%                 case 0.1
%                     r_b_0 = -0.007046;
%                     r_a_0 = 0.0281971;
%                 case 0.15
%                     converged = true;
%                     r_b_0 = -0.009953;
%                     r_a_0 = 0.034458;
%             end
%         case 5
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.038573;
%                     r_a_0 = 0.008789;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.039105;
%                     r_a_0 = 0.009146;
%                 case 0.02
%                     converged = true;
%                     r_b_0 = -0.040288;
%                     r_a_0 = 0.010153;
%                 case 0.05
%                     r_b_0 = -0.043872;
%                     r_a_0 = 0.020275;
%                 case 0.1
%                     r_b_0 = -0.051544;
%                     r_a_0 = 0.030421;
%                 case 0.15
%                     converged = true;
%                     r_b_0 = -0.061004;
%                     r_a_0 = 0.045998;
%             end
%         case 10
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.070741;
%                     r_a_0 = 0.005279;
%                 case 0.01
%                     converged = true;
%                     r_b_0 = -0.071900;
%                     r_a_0 = 0.005972;
%                 case 0.02
%                     r_b_0 = -0.074081;
%                     r_a_0 = 0.010808;
%                 case 0.05
%                     converged = true;
%                     r_b_0 = -0.082419;
%                     r_a_0 = 0.02236;
%                 case 0.1
%                     r_b_0 = -0.092600;
%                     r_a_0 = 0.039868;
%                 case 0.15
%                     converged = true;
%                     r_b_0 = -0.1099149;
%                     r_a_0 = 0.059595;
%             end
%         case 20
%             switch p.sigma_r
%                 case 0
%                     converged = true;
%                     r_b_0 = -0.092393;
%                     r_a_0 = 0.004989;
%                 case 0.01
%                     r_b_0 = -0.094769;
%                     r_a_0 = 0.005029;
%                 case 0.02
%                     r_b_0 = -0.100461;
%                     r_a_0 = 0.006972;
%                 case 0.05
%                     r_b_0 = -0.110461;
%                     r_a_0 = 0.011272;
%                 case 0.1
%                     r_b_0 = -0.114661;
%                     r_a_0 = 0.04572;
%                 case 0.15
%                     r_b_0 = -0.127261;
%                     r_a_0 = 0.075572;
%             end
%     end
% end
% 
% calibrator = solver.Calibrator(runopts, p, "r_b, r_a");
% x0 = calibrator.create_initial_condition([r_b_0, r_a_0]);

% opts = optimoptions('fsolve', 'MaxFunctionEvaluations', 400, 'MaxIterations', 600);
% fsolve(calibrator.objective, x0, opts);

% end

%% r_b, r_a calibration for rho calibrated to RA = 5
% calibrator = solver.Calibrator(runopts, p, "r_b, r_a");


% % r_a, rho calibration
% calibrator = solver.Calibrator(runopts, p, "r_a, rho");
% x0 = calibrator.create_initial_condition([0.066513, 0.142622]);
% fsolve(calibrator.objective, x0);

% % rho calibration
% calibrator = solver.Calibrator(runopts, p, "rho");
% x0 = calibrator.create_initial_condition(0.035);
% fsolve(calibrator.objective, x0);

% final run
p.set("NoRisk", 1);
p.set("ComputeMPCS", 1);
% p.set("ComputeMPCS_news", 1);
p.set("SaveResults", 1);
tic
stats = main(runopts, p);
toc
