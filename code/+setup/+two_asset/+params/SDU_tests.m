function outparams = SDU_tests(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,.

    %%--------------------------------------------------------------
    % BASELINE
    % --------------------------------------------------------------
    ii = 1;
    params(ii).name = sprintf('baseline (log utility)'); 
    params(ii).OneAsset = 0;
    params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
    params(ii).chi0 = 0;
    params(ii).chi1 = 0.15;
    params(ii).chi2 = 0.25;
    params(ii).a_lb = 0.25;
    params(ii).rhoL = 0.005;
    params(ii).rho = 0.015440584992491;
    params(ii).riskaver = 1;
    params(ii).invies = 1;
    params(ii).SDU = 0;
    params(ii).r_a = 0.0190643216;
    params(ii).delta_HJB = 1e4;
    params(ii).maxit_HJB = 1e6;
    params(ii).NoRisk = 0;
    params(ii).nb = 120;
    params(ii).nb_KFE = 120;

    %%--------------------------------------------------------------
    % WITH RETURNS RISK
    % --------------------------------------------------------------
    % illiquid_returns = linspace(0.03/4, 0.12/4, 10);
    risk_avers = [1, 2, 5, 10, 20];
    sdrs = [0, 0.01, 0.02, 0.05, 0.1, 0.15];

    ii = 2;
    for risk_aver = risk_avers
        for sd_r = sdrs
            params(ii).name = sprintf('SDU with riskaver%f, sigma_r%f', risk_aver, sd_r); 
            params(ii).OneAsset = 0;
            params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
            params(ii).chi0 = 0;
            params(ii).chi1 = 0.15;
            params(ii).chi2 = 0.25;
            params(ii).a_lb = 0.25;
            params(ii).rho = 0.015440584992491;
            params(ii).riskaver = risk_aver;
            params(ii).invies = 1;
            params(ii).SDU = 1;
            params(ii).maxit_HJB = 1e6;
            params(ii).sigma_r = sd_r;
            params(ii).retrisk_KFE = 0;
            params(ii).NoRisk = 0;
            params(ii).delta_HJB = 10;
            params(ii).delta_KFE = 10;
            params(ii).nb = 120;
            params(ii).nb_KFE = 120;
            
            if risk_aver == 1
                params(ii).SDU = 0;
            end

            if risk_aver == 10
                params(ii).delta_KFE = 1;
                params(ii).delta_HJB = 1;
            elseif risk_aver == 20
                params(ii).delta_KFE = 0.5;
                params(ii).delta_HJB = 0.1;
            end

            ii = ii + 1;
        end
    end

    %% DO NOT CHANGE BELOW

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = setup.two_asset.ParamsTwoAsset(runopts,chosen_param);

end