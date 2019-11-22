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

runopts.Server = 1; % sets IterateRho=1,fast=0,param_index=slurm env var
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


%% r_b, r_a calibration for rho calibrated to RA = 1, original adj costs

if p.chi1 == 0.15
    
if p.invies == 1
    switch p.riskaver
        case 1
            switch p.sigma_r
                case 0
                    r_b_0 = 0.006;
                    r_a_0 = 0.022866;
                case 0.01
                    r_b_0 = 0.005006;
                    r_a_0 = 0.022877;
                case 0.02
                    r_b_0 = 0.005026;
                    r_a_0 = 0.022909;
                case 0.05
                    r_b_0 = 0.005154;
                    r_a_0 = 0.02313;
                case 0.1
                    r_b_0 = 0.005539;
                    r_a_0 = 0.023891;
                case 0.15
                    r_b_0 = 0.006049;
                    r_a_0 = 0.025091;
            end
        case 2
            switch p.sigma_r
                case 0
                    r_b_0 = -0.009868;
                    r_a_0 = 0.01987;
                case 0.01
                    r_b_0 = -0.009856;
                    r_a_0 = 0.022877;
                case 0.02
                    r_b_0 = -0.00982;
                    r_a_0 = 0.020134;
                case 0.05
                    r_b_0 = -0.00975;
                    r_a_0 = 0.021134;
                case 0.1
                    r_b_0 = -0.009709;
                    r_a_0 = 0.024941;
                case 0.15
                    r_b_0 = -0.00961;
                    r_a_0 = 0.028941;
            end
        case 5
            switch p.sigma_r
                case 0
                    r_b_0 = -0.051442;
                    r_a_0 = 0.009666;
                case 0.01
                    r_b_0 = -0.051958;
                    r_a_0 = 0.009898;
                case 0.02
                    r_b_0 = -0.05323;
                    r_a_0 = 0.010561;
                case 0.05
                    r_b_0 = -0.057679;
                    r_a_0 = 0.014174;
                case 0.1
                    r_b_0 = -0.065535;
                    r_a_0 = 0.022612;
                case 0.15
                    r_b_0 = -0.07403;
                    r_a_0 = 0.03228;
            end
        case 10
            switch p.sigma_r
                case 0
                    r_b_0 = -0.092879;
                    r_a_0 = 0.00526;
                case 0.01
                    r_b_0 = -0.094133;
                    r_a_0 = 0.005732;
                case 0.02
                    r_b_0 = -0.096433;
                    r_a_0 = 0.006332;
                case 0.05
                    r_b_0 = -0.099433;
                    r_a_0 = 0.006832;
                case 0.1
                    r_b_0 = -0.105433;
                    r_a_0 = 0.007332;
                case 0.15
                    r_b_0 = -0.111433;
                    r_a_0 = 0.008032;
            end
        case 20
            switch p.sigma_r
                case 0
                    r_b_0 = -0.124384;
                    r_a_0 = 0.004371;
                case 0.01
                    r_b_0 = -0.128384;
                    r_a_0 = 0.006371;
                case 0.02
                    r_b_0 = -0.133384;
                    r_a_0 = 0.007371;
                case 0.05
                    r_b_0 = -0.144384;
                    r_a_0 = 0.008371;
                case 0.1
                    r_b_0 = -0.152;
                    r_a_0 = 0.009471;
                case 0.15
                    r_b_0 = -0.1689;
                    r_a_0 = 0.012;
            end
    end
else
    switch p.riskaver
        case 1.01
            switch p.sigma_r
                case 0
                    r_b_0 = 0.016998;
                    r_a_0 = 0.027251;
                case 0.01
                    r_b_0 = 0.017017;
                    r_a_0 = 0.027289;
                case 0.02
                    r_b_0 = 0.017073;
                    r_a_0 = 0.027399;
                case 0.05
                    r_b_0 = 0.017508;
                    r_a_0 = 0.028097;
                case 0.1
                    r_b_0 = 0.018308;
                    r_a_0 = 0.029097;
                case 0.15
                    r_b_0 = 0.019308;
                    r_a_0 = 0.030597;
            end
        case 2
            switch p.sigma_r
                case 0
                    r_b_0 = 0.00373;
                    r_a_0 = 0.024604;
                case 0.01
                    r_b_0 = 0.003824;
                    r_a_0 = 0.024692;
                case 0.02
                    r_b_0 = 0.004089;
                    r_a_0 = 0.024949;
                case 0.05
                    r_b_0 = 0.005161;
                    r_a_0 = 0.026571;
                case 0.1
                    r_b_0 = 0.006461;
                    r_a_0 = 0.0281971;
                case 0.15
                    r_b_0 = 0.007861;
                    r_a_0 = 0.0321971;
            end
        case 5
            switch p.sigma_r
                case 0
                    r_b_0 = -0.036748;
                    r_a_0 = 0.014726;
                case 0.01
                    r_b_0 = -0.036771;
                    r_a_0 = 0.014992;
                case 0.02
                    r_b_0 = -0.036853;
                    r_a_0 = 0.015755;
                case 0.05
                    r_b_0 = -0.038872;
                    r_a_0 = 0.020275;
                case 0.1
                    r_b_0 = -0.047149;
                    r_a_0 = 0.031595;
                case 0.15
                    r_b_0 = -0.059149;
                    r_a_0 = 0.043595;
            end
        case 10
            switch p.sigma_r
                case 0
                    r_b_0 = -0.077692;
                    r_a_0 = 0.010993;
                case 0.01
                    r_b_0 = -0.077745;
                    r_a_0 = 0.011551;
                case 0.02
                    r_b_0 = -0.07771;
                    r_a_0 = 0.013081;
                case 0.05
                    r_b_0 = -0.082419;
                    r_a_0 = 0.02136;
                case 0.1
                    r_b_0 = -0.092419;
                    r_a_0 = 0.031595;
                case 0.15
                    r_b_0 = -0.109149;
                    r_a_0 = 0.043595;
            end
        case 20
            switch p.sigma_r
                case 0
                    r_b_0 = -0.1018;
                    r_a_0 = 0.010889;
                case 0.01
                    r_b_0 = -0.104769;
                    r_a_0 = 0.012029;
                case 0.02
                    r_b_0 = -0.108461;
                    r_a_0 = 0.014972;
                case 0.05
                    r_b_0 = -0.110461;
                    r_a_0 = 0.016272;
                case 0.1
                    r_b_0 = -0.114661;
                    r_a_0 = 0.019372;
                case 0.15
                    r_b_0 = -0.127261;
                    r_a_0 = 0.023572;
            end
    end
end

calibrator = solver.Calibrator(runopts, p, "r_b, r_a");
x0 = calibrator.create_initial_condition([r_b_0, r_a_0]);

opts = optimoptions('fsolve', 'MaxFunctionEvaluations', 400, 'MaxIterations', 600);
fsolve(calibrator.objective, x0, opts);

end

%% r_b, r_a calibration for rho calibrated to RA = 1, NEW adj costs

if p.chi1 ~= 0.15
    
if p.invies == 1
    switch p.riskaver
        case 1
            switch p.sigma_r
                case 0
                    r_b_0 = 0.005;
                    r_a_0 = 0.01979;
                case 0.01
                    r_b_0 = 0.005006;
                    r_a_0 = 0.020877;
                case 0.02
                    r_b_0 = 0.005026;
                    r_a_0 = 0.020909;
                case 0.05
                    r_b_0 = 0.005154;
                    r_a_0 = 0.02113;
                case 0.1
                    r_b_0 = 0.005539;
                    r_a_0 = 0.021891;
                case 0.15
                    r_b_0 = 0.006049;
                    r_a_0 = 0.022091;
            end
        case 2
            switch p.sigma_r
                case 0
                    r_b_0 = -0.009868;
                    r_a_0 = 0.01987;
                case 0.01
                    r_b_0 = -0.009856;
                    r_a_0 = 0.022877;
                case 0.02
                    r_b_0 = -0.00982;
                    r_a_0 = 0.020134;
                case 0.05
                    r_b_0 = -0.00975;
                    r_a_0 = 0.021134;
                case 0.1
                    r_b_0 = -0.009709;
                    r_a_0 = 0.024941;
                case 0.15
                    r_b_0 = -0.00961;
                    r_a_0 = 0.028941;
            end
        case 5
            switch p.sigma_r
                case 0
                    r_b_0 = -0.051442;
                    r_a_0 = 0.009666;
                case 0.01
                    r_b_0 = -0.051958;
                    r_a_0 = 0.009898;
                case 0.02
                    r_b_0 = -0.05323;
                    r_a_0 = 0.010561;
                case 0.05
                    r_b_0 = -0.057679;
                    r_a_0 = 0.014174;
                case 0.1
                    r_b_0 = -0.065535;
                    r_a_0 = 0.022612;
                case 0.15
                    r_b_0 = -0.07403;
                    r_a_0 = 0.03228;
            end
        case 10
            switch p.sigma_r
                case 0
                    r_b_0 = -0.092879;
                    r_a_0 = 0.00526;
                case 0.01
                    r_b_0 = -0.094133;
                    r_a_0 = 0.005732;
                case 0.02
                    r_b_0 = -0.096433;
                    r_a_0 = 0.006332;
                case 0.05
                    r_b_0 = -0.099433;
                    r_a_0 = 0.006832;
                case 0.1
                    r_b_0 = -0.105433;
                    r_a_0 = 0.007332;
                case 0.15
                    r_b_0 = -0.111433;
                    r_a_0 = 0.008032;
            end
        case 20
            switch p.sigma_r
                case 0
                    r_b_0 = -0.124384;
                    r_a_0 = 0.004371;
                case 0.01
                    r_b_0 = -0.128384;
                    r_a_0 = 0.006371;
                case 0.02
                    r_b_0 = -0.133384;
                    r_a_0 = 0.007371;
                case 0.05
                    r_b_0 = -0.144384;
                    r_a_0 = 0.008371;
                case 0.1
                    r_b_0 = -0.152;
                    r_a_0 = 0.009471;
                case 0.15
                    r_b_0 = -0.1689;
                    r_a_0 = 0.012;
            end
    end
else
    switch p.riskaver
        case 1.01
            switch p.sigma_r
                case 0
                    r_b_0 = 0.016998;
                    r_a_0 = 0.027251;
                case 0.01
                    r_b_0 = 0.017017;
                    r_a_0 = 0.027289;
                case 0.02
                    r_b_0 = 0.017073;
                    r_a_0 = 0.027399;
                case 0.05
                    r_b_0 = 0.017508;
                    r_a_0 = 0.028097;
                case 0.1
                    r_b_0 = 0.018308;
                    r_a_0 = 0.029097;
                case 0.15
                    r_b_0 = 0.019308;
                    r_a_0 = 0.030597;
            end
        case 2
            switch p.sigma_r
                case 0
                    r_b_0 = 0.00373;
                    r_a_0 = 0.024604;
                case 0.01
                    r_b_0 = 0.003824;
                    r_a_0 = 0.024692;
                case 0.02
                    r_b_0 = 0.004089;
                    r_a_0 = 0.024949;
                case 0.05
                    r_b_0 = 0.005161;
                    r_a_0 = 0.026571;
                case 0.1
                    r_b_0 = 0.006461;
                    r_a_0 = 0.0281971;
                case 0.15
                    r_b_0 = 0.007861;
                    r_a_0 = 0.0321971;
            end
        case 5
            switch p.sigma_r
                case 0
                    r_b_0 = -0.036748;
                    r_a_0 = 0.014726;
                case 0.01
                    r_b_0 = -0.036771;
                    r_a_0 = 0.014992;
                case 0.02
                    r_b_0 = -0.036853;
                    r_a_0 = 0.015755;
                case 0.05
                    r_b_0 = -0.038872;
                    r_a_0 = 0.020275;
                case 0.1
                    r_b_0 = -0.047149;
                    r_a_0 = 0.031595;
                case 0.15
                    r_b_0 = -0.059149;
                    r_a_0 = 0.043595;
            end
        case 10
            switch p.sigma_r
                case 0
                    r_b_0 = -0.077692;
                    r_a_0 = 0.010993;
                case 0.01
                    r_b_0 = -0.077745;
                    r_a_0 = 0.011551;
                case 0.02
                    r_b_0 = -0.07771;
                    r_a_0 = 0.013081;
                case 0.05
                    r_b_0 = -0.082419;
                    r_a_0 = 0.02136;
                case 0.1
                    r_b_0 = -0.092419;
                    r_a_0 = 0.031595;
                case 0.15
                    r_b_0 = -0.109149;
                    r_a_0 = 0.043595;
            end
        case 20
            switch p.sigma_r
                case 0
                    r_b_0 = -0.1018;
                    r_a_0 = 0.010889;
                case 0.01
                    r_b_0 = -0.104769;
                    r_a_0 = 0.012029;
                case 0.02
                    r_b_0 = -0.108461;
                    r_a_0 = 0.014972;
                case 0.05
                    r_b_0 = -0.110461;
                    r_a_0 = 0.016272;
                case 0.1
                    r_b_0 = -0.114661;
                    r_a_0 = 0.019372;
                case 0.15
                    r_b_0 = -0.127261;
                    r_a_0 = 0.023572;
            end
    end
end

calibrator = solver.Calibrator(runopts, p, "r_b, r_a");
x0 = calibrator.create_initial_condition([r_b_0, r_a_0]);

opts = optimoptions('fsolve', 'MaxFunctionEvaluations', 400, 'MaxIterations', 600);
fsolve(calibrator.objective, x0, opts);

end

%% r_b, r_a calibration for rho calibrated to RA = 5
% calibrator = solver.Calibrator(runopts, p, "r_b, r_a");
% if p.invies == 1
%     switch p.riskaver
%         case 1
%             x0 = calibrator.create_initial_condition([0.025, 0.06]);
%         case 2
%             x0 = calibrator.create_initial_condition([-0.04, 0.0001]);
%         case 5
%             x0 = calibrator.create_initial_condition([-0.06, -0.001]);
%         case 10
%             x0 = calibrator.create_initial_condition([-0.1, -0.002]);
%         case 20
%             x0 = calibrator.create_initial_condition([-0.12, -0.005]);
%     end
% else
%     switch p.riskaver
%         case 1.01
%             x0 = calibrator.create_initial_condition([0.005, 0.014]);
%         case 2
%             x0 = calibrator.create_initial_condition([-0.035, 0.005]);
%         case 5
%             x0 = calibrator.create_initial_condition([-0.09, -0.0005]);
%         case 10
%             x0 = calibrator.create_initial_condition([-0.125, -0.006]);
%         case 20
%             x0 = calibrator.create_initial_condition([-0.15, -0.01]);
%     end
% end
% fsolve(calibrator.objective, x0);

%% r_a, rho calibration
% calibrator = solver.Calibrator(runopts, p, "r_a, rho");
% x0 = calibrator.create_initial_condition([0.043541, 0.0319]);
% fsolve(calibrator.objective, x0);

%% rho calibration
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
