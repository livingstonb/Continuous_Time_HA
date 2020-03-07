function outparams = params_one_asset(param_opts, param_index)

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
    params{ii}.n_mpcsim = 5e5;
    params{ii}.nb = 200;
    params{ii}.nb_KFE = 200;
    ii = ii + 1;

    %% --------------------------------------------------------------------
    % HOUSEKEEPING, DO NOT CHANGE
    % ---------------------------------------------------------------------

    % Use runopts.param_index to choose which specification to select
    chosen_param = params{param_index};

    % Create Params object
    outparams = HACTLib.model_objects.Params(chosen_param);
end