function outparams = params_one_asset(param_opts, param_index)

    dollars = [-1, -500, -5000, 1, 500, 5000];
    shocks = dollars ./ 72000;

	%% --------------------------------------------------------------------
    % CREATE PARAMETERIZATIONS HERE
    % ---------------------------------------------------------------------
    ii = 1;
    params{ii} = param_opts;
	params{ii}.name = 'calib_to_median_wealth'; 
    params{ii}.OneAsset = 1;
    params{ii}.income_dir = 'continuous_a';
    params{ii}.Bequests = 1;
    params{ii}.rho = 0.003902728727572;
    params{ii}.n_mpcsim = 5e5;
    params{ii}.nb = 200;
    params{ii}.nb_KFE = 200;
    ii = ii + 1;

    params{ii} = param_opts;
    params{ii}.name = 'baseline'; 
    params{ii}.OneAsset = 1;
    params{ii}.income_dir = 'continuous_a';
    params{ii}.Bequests = 1;
    params{ii}.rho = 0.003902728727572;
    params{ii}.nb = 250;
    params{ii}.nb_KFE = 250;
    params{ii}.mpc_shocks = shocks;
    params{ii}.b_gcurv_pos = 0.3;
    params{ii}.mpc_shocks_dollars = [];
    params{ii}.calibration_vars = {'rho'};
    params{ii}.calibration_stats = {'totw'};
    params{ii}.calibration_targets = 3.5;
    params{ii}.calibration_bounds = {[0.003, 0.0042]};
    params{ii}.bmax = 50;
    ii = ii + 1;
    
    params{ii} = param_opts;
    params{ii}.name = 'cont_b'; 
    params{ii}.OneAsset = 1;
    params{ii}.income_dir = 'continuous_b';
    params{ii}.Bequests = 1;
    params{ii}.rho = 0.003902728727572;
    params{ii}.nb = 250;
    params{ii}.nb_KFE = 250;
    params{ii}.mpc_shocks = shocks;
    params{ii}.b_gcurv_pos = 0.3;
    params{ii}.mpc_shocks_dollars = [];
    params{ii}.calibration_vars = {'rho'};
    params{ii}.calibration_stats = {'totw'};
    params{ii}.calibration_targets = 3.5;
    params{ii}.calibration_bounds = {[0.0042, 0.005]};
    params{ii}.bmax = 50;
    ii = ii + 1;

    
    %% --------------------------------------------------------------------
    % HOUSEKEEPING, DO NOT CHANGE
    % ---------------------------------------------------------------------

    % Use runopts.param_index to choose which specification to select
    chosen_param = params{param_index};

    % Create Params object
    outparams = HACTLib.model_objects.Params(chosen_param);
end