function [outparams, n] = params_one_asset(param_opts)
    import setup.params.scf2019struct
    scf = scf2019struct();

	shocks = [-1, -500, -5000, 1, 500, 5000];

    shared_params = param_opts;
    shared_params.mpc_shocks = shocks / (scf.quarterly_earnings * 4);
    shared_params.numeraire_in_dollars = (scf.quarterly_earnings * 4);
%     shared_params.nb = 250;
%     shared_params.nb_KFE = 250;
    shared_params.nb = 150;
    shared_params.nb_KFE = 150;

    shared_params.bgrid_term1_weight = 0;
    shared_params.bgrid_term1_curv = 1;
    shared_params.b_gcurv_pos = 0.2;
    shared_params.OneAsset = true;
    shared_params.Bequests = true;
    shared_params.no_transitory_incrisk = false;

    shared_params.bmax = 500;
    shared_params.rho = 0.005;

    shared_params.calibration_vars = {'rho'};
    shared_params.calibration_stats = {'totw'};
    shared_params.calibration_targets = 3.5;
    shared_params.calibration_bounds = {[0.003, 0.01]};
    
    median_calibration = shared_params;
    median_calibration.calibration_vars = {'rho'};
    median_calibration.calibration_bounds = {[0.001, 0.05]};
    median_calibration.calibration_backup_x0 = {[0.004]};
    median_calibration.calibration_stats = {'median_totw'};
    median_calibration.calibration_targets = [scf.median_totw];
    median_calibration.calibration_scales = [1];

	%% --------------------------------------------------------------------
    % CREATE PARAMETERIZATIONS HERE
    % ---------------------------------------------------------------------
    ii = 1;
%     params{ii} = shared_params;
% 	params{ii}.name = 'calib_to_median_wealth'; 
%     params{ii}.income_dir = 'continuous_a';
%     params{ii}.calibration_vars = {'rho'};
%     params{ii}.calibration_stats = {'median_totw'};
%     params{ii}.calibration_targets = 1.7;
%     params{ii}.calibration_bounds = {[0.002, 0.01]};
%     ii = ii + 1;

    params{ii} = median_calibration;
    params{ii}.name = 'baseline'; 
    params{ii}.income_dir = 'continuous_a';
    params{ii}.IncomeDescr = 'cont_a';
    params{ii}.rho = 0.003902728727572;
    ii = ii + 1;
    
    params{ii} = median_calibration;
    params{ii}.name = 'cont_b'; 
    params{ii}.income_dir = 'continuous_b';
    params{ii}.IncomeDescr = 'cont_b';
    params{ii}.rho = 0.003902728727572;
    ii = ii + 1;
    
    params{ii} = median_calibration;
    params{ii}.name = 'cont_a_no_meas_error'; 
    params{ii}.income_dir = 'continuous_a/no_measurement_error';
    params{ii}.IncomeDescr = 'cont_a_no_meas_error';
    params{ii}.rho = 0.003902728727572;
    ii = ii + 1;

    params{ii} = median_calibration;
    params{ii}.name = 'cont_a_20pc_meas_error'; 
    params{ii}.income_dir = 'continuous_a/measurement_error_20pc';
    params{ii}.IncomeDescr = 'cont_a_20pc_meas_error';
    params{ii}.rho = 0.003902728727572;
    ii = ii + 1;
    
    params{ii} = median_calibration;
    params{ii}.name = 'cont_a_33pc_meas_error'; 
    params{ii}.income_dir = 'continuous_a/measurement_error_33pc';
    params{ii}.IncomeDescr = 'cont_a_33pc_meas_error';
    params{ii}.rho = 0.003902728727572;
    ii = ii + 1;
    
    params{ii} = median_calibration;
    params{ii}.name = 'cont_a_50pc_meas_error'; 
    params{ii}.income_dir = 'continuous_a/measurement_error_50pc';
    params{ii}.IncomeDescr = 'cont_a_50pc_meas_error';
    params{ii}.rho = 0.003902728727572;
    ii = ii + 1;
    
    %% --------------------------------------------------------------------
    % HOUSEKEEPING, DO NOT CHANGE
    % ---------------------------------------------------------------------
    n = numel(params);
    outparams = params{param_opts.param_index};
end