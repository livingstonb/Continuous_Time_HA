function outparams = endog_labor_tests(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,...
    
    % Any parameters not set here will take on their default values as in
    % Classes/Params.m
    
    % An example for parameterization 1 is:
    % params(1).name = 'Parameterization 1';
    % params(1).DirIncomeProcess = 'IncomeGrids/continuous_b';
    % params(1).nb = 50;
    
    %% --------------------------------------------------------------------
    % CREATE PARAMETERIZATIONS HERE
    % ---------------------------------------------------------------------
    
    ii = 1;

    % one asset
    params(ii).name = 'one_asset_endog_labor';
    params(ii).OneAsset = 1;
    params(ii).income_dir = 'twopoint_3_5';
    params(ii).targetAY = 3.5;
    params(ii).Bequests = 0;
    params(ii).deathrate = 0;
    params(ii).endogenous_labor = true;
    params(ii).nb = 50;
    params(ii).nb_KFE = 50;
    params(ii).MPL = 1/3;
    params(ii).transfer = 0.0081 * 2;
    ii = ii + 1;

    % two asset
    params(ii).name = 'two_asset_endog_labor';
    params(ii).OneAsset = 0;
    params(ii).income_dir = 'twopoint_3_5';
    params(ii).Bequests = 0;
    params(ii).deathrate = 0;
    params(ii).chi0 = 0;
    params(ii).chi1 = 0.15;
    params(ii).chi2 = 0.25;
    params(ii).a_lb = 0.25;
    params(ii).rho = 0.01;
    params(ii).r_a = 0.0190643216;
    params(ii).nb = 50;
    params(ii).nb_KFE = 50;
    params(ii).na = 50;
    params(ii).na_KFE = 50;
    ii = ii + 1;


    %% --------------------------------------------------------------------
    % HOUSEKEEPING, DO NOT CHANGE
    % ---------------------------------------------------------------------
    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACTLib.model_objects.Params(runopts,chosen_param);

end
