%% ONE-ASSET HOUSEHOLD MODEL

% This is the main script for this code repository
% HA model with a liquid asset and costly consumption adjustment

% Prior to running this script:

% (1) Set options in the section below. Directory only needs to be set for
% either the server or the local directory, depending on whether
% runopts.Server = 0 or 1

% (2) Modify the parameters script 'get_params.m' and make sure that 
% runopts.mode is equal to 'get_params'. Note that all parameter defaults
% are set in the class file Model_Setup/Params.m, and get_params.m overrides
% these defaults. Any parameters not in get_params.m are set to their
% defaults. See the attributes of Model_Setup/Params.m for a list of all
% parameters.

% (3) Set runopts.param_index equal to the index of the parameterization
% you would like to run in the parameters file (1,2,...).

% RUNNING ON THE SERVER: To run in batch on the server, use 
% batch/server.sbatch as a template. That script sends an array to SLURM 
% that runs all of the requested parameters in get_params.m. Output files
% are stored in the Output directory

clear all
warning('off','MATLAB:nearlySingularMatrix')

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------

runopts.Server      = 0; % sets param_index=slurm env var
runopts.IterateRho  = 0;
runopts.fast = 1; % use small grid for  debugging
runopts.ComputeMPCS = 1;
runopts.ComputeMPCS_news = 1;
runopts.SimulateMPCS = 1;
runopts.SimulateMPCS_news = 1;

% Select which parameterization to run
% (ignored when runops.Server = 1)
runopts.param_index = 1;

% Location of Continuous_Two_Asset directory
runopts.serverdir = '/home/livingstonb/GitHub/Continuous_Time_HA/';
runopts.localdir = '/Users/brianlivingston/Documents/GitHub/Continuous_Time_HA/';

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE
% -------------------------------------------------------------------------

if runopts.Server == 0
	runopts.direc = runopts.localdir;
else
	runopts.direc = runopts.serverdir;
	runopts.param_index = str2num(getenv('SLURM_ARRAY_TASK_ID'));
	
	runopts.fast = 0;
	runopts.IterateRho = 1;
end
    
% for saving
runopts.suffix = num2str(runopts.param_index);

% directory to save output
runopts.savedir = [runopts.direc 'output/con_effort/'];
% temp directory
runopts.temp = [runopts.direc 'temp/con_effort/run' runopts.suffix '/'];

addpath([runopts.direc 'code']);
addpath([runopts.direc 'code/factorization_lib']);

mkdir(runopts.temp);
mkdir(runopts.savedir);
addpath(runopts.temp);
addpath(runopts.savedir);


% check that specified directories exist
if ~exist(runopts.direc,'dir')
    error([runopts.direc ' does not exist'])
end

cd(runopts.direc)

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION FILE
% -------------------------------------------------------------------------

tic
[stats,p,grdKFE,KFE] = main_con_effort(runopts);
toc

delete([runopts.temp '*'])
rmdir(runopts.temp)
