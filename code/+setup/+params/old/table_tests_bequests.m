function outparams = table_tests_bequests(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,.
    
    
    %% adjusing (A/Y,B/Y), fixed adj costs
    raBASELINE = 0.015132;
    
    % ras = [raBASELINE,0.0129226573755814,0.0304148974027774];
    ras = raBASELINE;
        
    % calibrate r_a to match B/Y in each case
    ii = 1;
    for ra = ras
        params(ii).name = sprintf('contB, fixed_adj_costs, r_a=%f',ra); 
        params(ii).OneAsset = 0;
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).chi0 = 0;
        params(ii).chi1 = 0.15;
        params(ii).chi2 = 0.25;
        params(ii).a_lb = 0.25;
        params(ii).rhoL = 0.02;
        params(ii).rho = 0.015440584992491;
        params(ii).Bequests = 1;
        params(ii).maxit_KFE = 2e4;
        params(ii).delta_KFE = 1e3;
        params(ii).r_a = ra;
        
        ii = ii + 1;
    end
    
    %% varying chi1 and chi2
    ii = 101;
    for chi1 = [0.1,0.4]
        params(ii).name = sprintf('chi1=%f',chi1);
        params(ii).OneAsset = 0;
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).chi0 = 0;
        params(ii).chi1 = chi1;
        params(ii).chi2 = 0.25;
        params(ii).a_lb = 0.25;
        params(ii).rhoL = 0.02;
        params(ii).maxit_KFE = 2e4;
        params(ii).delta_KFE = 1e3;
        params(ii).r_a = raBASELINE;
        
        ii = ii + 1;
    end
    
    for chi2 = [0.1,0.15,0.5,1]
        params(ii).name = sprintf('chi1=%f',chi1);
        params(ii).OneAsset = 0;
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).chi0 = 0;
        params(ii).chi1 = 0.15;
        params(ii).chi2 = chi2;
        params(ii).a_lb = 0.25;
        params(ii).rhoL = 0.02;
        params(ii).maxit_KFE = 2e4;
        params(ii).delta_KFE = 1e3;
        params(ii).r_a = raBASELINE;
        
        ii = ii + 1;
    end
    
    %% varying r_a and r_b
    
    ras_annual = 0.05:0.01:0.1;
    for i = 1:6
        ras_quarterly(i) = (1+ras_annual(i)).^(1/4) - 1;
    end
    
    ii = 201;
    for ra = ras_quarterly
        params(ii).name = sprintf('annual r_a=%f',(1+ra)^4-1);
        params(ii).OneAsset = 0;
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).chi0 = 0;
        params(ii).chi1 = 0.15;
        params(ii).chi2 = 0.25;
        params(ii).a_lb = 0.25;
        params(ii).rhoL = 0.02;
        params(ii).maxit_KFE = 2e4;
        params(ii).delta_KFE = 1e3;
        params(ii).r_a = ra;
        
        ii = ii + 1;
    end
    
    rbs_annual = [0,0.05];
    for i = 1:2
        rbs_quarterly(i) = (1+rbs_annual(i)).^(1/4) - 1;
    end
    
    for rb = rbs_quarterly
        params(ii).name = sprintf('annual r_b=%f',(1+rb)^4-1);
        params(ii).OneAsset = 0;
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).chi0 = 0;
        params(ii).chi1 = 0.15;
        params(ii).chi2 = 0.25;
        params(ii).a_lb = 0.25;
        params(ii).rhoL = 0.02;
        params(ii).r_b = rb;
        params(ii).r_a = 0.0190643216;
        params(ii).maxit_KFE = 2e4;
        params(ii).delta_KFE = 1e3;
        
        ii = ii + 1;
    end
    
    %% recalibrate (rho,r_a) to match B/Y=0.5 with cont_a
    ii = 301;
    ras = 0.0173096090154433;
    for ra = ras
        params(ii).name = 'cont_a, recalibrated'; 
        params(ii).OneAsset = 0;
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_a';
        params(ii).chi0 = 0;
        params(ii).chi1 = 0.15;
        params(ii).chi2 = 0.25;
        params(ii).a_lb = 0.25;
        params(ii).rhoL = 0.008;
        params(ii).r_a = ra;
        params(ii).maxit_KFE = 2e4;
        params(ii).delta_KFE = 1e3;
        
        ii = ii + 1;
    end
    
    %% run cont_a without re-calibrating (rho,r_a)
    ii = 401;
    params(ii).name = sprintf('conta, no recalibration'); 
    params(ii).OneAsset = 0;
    params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_a';
    params(ii).chi0 = 0;
    params(ii).chi1 = 0.15;
    params(ii).chi2 = 0.25;
    params(ii).a_lb = 0.25;
    params(ii).rhoL = 0.02;
    params(ii).r_a = 0.01906432159025;
    params(ii).rho = 0.015440584980253;
    params(ii).maxit_KFE = 2e4;
    params(ii).delta_KFE = 1e3;

    %% risk aversion coeff tests
    ii = 501;

    for riskaver = [0.5,2,4,6]
        params(ii).name = sprintf('risk aversion = %f',riskaver); 
        params(ii).OneAsset = 0;
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).chi0 = 0;
        params(ii).chi1 = 0.15;
        params(ii).chi2 = 0.25;
        params(ii).a_lb = 0.25;
        params(ii).rho = 0.02;
        params(ii).r_a = raBASELINE;
        params(ii).riskaver = riskaver;
        params(ii).maxit_KFE = 2e4;
        params(ii).delta_KFE = 1e3;
        
        if riskaver == 0.5
            params(ii).rhoL = 0.01;
        elseif riskaver == 2
            params(ii).rhoL = 0.017;
        elseif riskaver == 4
            params(ii).delta_KFE = 1e3;
            params(ii).rhoL = 0.024;
        elseif riskaver == 6
            params(ii).delta_KFE = 1e1;
            params(ii).rhoL = 0.04;
        end
        ii = ii + 1;
    end

    %% rho heterogeneity
    ii = 601;

    for w = [0.001,0.005]
        params(ii).name = sprintf('rho het, width = %f',w); 
        params(ii).OneAsset = 0;
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).chi0 = 0;
        params(ii).chi1 = 0.15;
        params(ii).chi2 = 0.25;
        params(ii).a_lb = 0.25;
        params(ii).rho = 0.005;
        params(ii).rhoL = 0.02;
        params(ii).r_a = raBASELINE;
        params(ii).rho_grid = [-w,0,w];
        params(ii).maxit_KFE = 2e4;
        params(ii).delta_KFE = 1e3;
        ii = ii + 1;
    end

    %% DO NOT CHANGE BELOW

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACTLib.model_objects.Params(runopts,chosen_param);

end