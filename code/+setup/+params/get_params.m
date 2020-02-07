function outparams = get_params(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,...
    
    % Any parameters not set here will take on their default values as in
    % HACTLib/+model_objects/ParamsDefaults.m
    
    % An example for parameterization 1 is:
    % params(1).name = 'Parameterization 1';
    % params(1).DirIncomeProcess = 'IncomeGrids/continuous_b';
    % params(1).nb = 50;
    
    %% --------------------------------------------------------------------
    % CREATE PARAMETERIZATIONS HERE
    % ---------------------------------------------------------------------
    
    i = 1;
    
    for target = [0.25]
        params(i).name = 'baseline_cont_a'; 
        params(i).OneAsset = 1;
        params(i).income_dir = 'continuous_a';
        params(i).targetAY = target;
        params(i).Bequests = 1;
        params(i).rho = 0.008;
        params(i).n_mpcsim = 5e5;
        i = i + 1;

        params(i).name = 'continuous_b';
        params(i).OneAsset = 1;
        params(i).income_dir = 'continuous_b';
        params(i).targetAY = target;
        params(i).Bequests = 1;
        params(i).n_mpcsim = 1e5;
        i = i + 1;
    end

    % two asset
    i = 101;
    params(i).name = 'two_asset';
    params(i).OneAsset = 0;
    params(i).income_dir = 'continuous_b';
    params(i).Bequests = 0;
    params(i).chi0 = 0;
    params(i).chi1 = 0.15;
    params(i).chi2 = 0.25;
    params(i).a_lb = 0.25;
    params(i).rhoL = 0.01;
    params(i).rho = 0.005;
    params(i).r_a = 0.0190643216;
    params(i).min_grid_spacing = 0.005;
    
    params(i).nb = 60;
    params(i).nb_pos = 60;
    params(i).nb_KFE = 60;
    params(i).nb_pos_KFE = 60;
    params(i).bmin = 0;

    %% discount rate heterogeneity test
    i = 201;
    params(i).name = 'two_asset_rho_heterogeneity';
    params(i).OneAsset = 0;
    params(i).income_dir = 'continuous_b';
    params(i).Bequests = 0;
    params(i).chi0 = 0;
    params(i).chi1 = 0.15;
    params(i).chi2 = 0.25;
    params(i).a_lb = 0.25;
    params(i).rhoL = 0.01;
    params(i).rho = 0.02;
    params(i).rho_grid = [-0.0005,0.0005]; % distance between them will be fixed
    params(i).r_a = 0.0190643216;

    params(i).nb = 30;
    params(i).nb_pos = 30;
    params(i).nb_KFE = 30;
    params(i).nb_pos_KFE = 30;
    
    i = 202;
    params(i).name = 'two_asset_no_rho_heterogeneity';
    params(i).OneAsset = 0;
    params(i).income_dir = 'continuous_b';
    params(i).Bequests = 0;
    params(i).chi0 = 0;
    params(i).chi1 = 0.15;
    params(i).chi2 = 0.25;
    params(i).a_lb = 0.25;
    params(i).rhoL = 0.01;
    params(i).rho = 0.008;
    params(i).rho_grid = 0;
    params(i).r_a = 0.0190643216;

    params(i).nb = 45;
    params(i).nb_pos = 45;
    params(i).nb_KFE = 45;
    params(i).nb_pos_KFE = 45;

    %% baseline one-asset with loose borrowing constraint
    i = 301;
    params(i).name = 'baseline_no_bc';
    params(i).OneAsset = 1;
    params(i).income_dir = 'continuous_b';
    params(i).targetAY = 3.5;
    params(i).Bequests = 1;
    params(i).n_mpcsim = 5e5;
    params(i).nb_pos = 500;
    params(i).nb_pos_KFE = 400;
    params(i).nb = 600;
    params(i).nb_KFE = 500;
    params(i).bmin = -1e10;
    i = i + 1;

    %% --------------------------------------------------------------------
    % HOUSEKEEPING, DO NOT CHANGE
    % ---------------------------------------------------------------------

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACTLib.model_objects.Params(runopts,chosen_param);

end
