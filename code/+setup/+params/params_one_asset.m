function [outparams, n] = params_one_asset(param_opts)

    shocks = [-1, -500, -5000, 1, 500, 5000];
    dollars = 72000;

    shared_params = param_opts;
    shared_params.mpc_shocks = shocks / dollars;
    shared_params.numeraire_in_dollars = dollars;
    shared_params.nb = 250;
    shared_params.nb_KFE = 250;

    shared_params.bgrid_term1_weight = 0.01;
    shared_params.bgrid_term1_curv = 0.8;
    shared_params.b_gcurv_pos = 0.1;
    shared_params.OneAsset = true;
    shared_params.Bequests = true;

    shared_params.bmax = 50;
    shared_params.rho = 0.005;

    shared_params.calibration_vars = {'rho'};
    shared_params.calibration_stats = {'totw'};
    shared_params.calibration_targets = 3.5;
    shared_params.calibration_bounds = {[0.003, 0.01]};

	%% --------------------------------------------------------------------
    % CREATE PARAMETERIZATIONS HERE
    % ---------------------------------------------------------------------
    ii = 1;
    params{ii} = shared_params;
	params{ii}.name = 'calib_to_median_wealth'; 
    params{ii}.income_dir = 'continuous_a';
    params{ii}.calibration_vars = {'rho'};
    params{ii}.calibration_stats = {'median_totw'};
    params{ii}.calibration_targets = 1.7;
    params{ii}.calibration_bounds = {[0.002, 0.01]};
    ii = ii + 1;

    params{ii} = shared_params;
    params{ii}.name = 'baseline'; 
    params{ii}.income_dir = 'continuous_a';
    params{ii}.rho = 0.003902728727572;
    ii = ii + 1;
    
    params{ii} = shared_params;
    params{ii}.name = 'cont_b'; 
    params{ii}.income_dir = 'continuous_b';
    params{ii}.rho = 0.003902728727572;
    ii = ii + 1;

    
    %% --------------------------------------------------------------------
    % HOUSEKEEPING, DO NOT CHANGE
    % ---------------------------------------------------------------------
    n = numel(params);
    outparams = params{param_opts.param_index};
end