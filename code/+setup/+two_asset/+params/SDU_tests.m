function outparams = SDU_tests(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,.

    ii = 1;
    params(ii).name = sprintf('SDU_test'); 
    params(ii).OneAsset = 0;
    params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
    params(ii).chi0 = 0;
    params(ii).chi1 = 0.15;
    params(ii).chi2 = 0.25;
    params(ii).a_lb = 0.25;
    params(ii).rhoL = 0.005;
    params(ii).rho = 0.05;
    params(ii).r_a = 0.06/4;
    params(ii).riskaver = 4;
    params(ii).invies = 1.2;
    params(ii).SDU = 1;
    params(ii).deathrate = 0;
    params(ii).delta_HJB = 1e4;
    params(ii).nb = 75;
    params(ii).nb_KFE = 75;
    params(ii).na = 75;
    params(ii).na_KFE = 75;

    %% DO NOT CHANGE BELOW

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = setup.two_asset.ParamsTwoAsset(runopts,chosen_param);

end