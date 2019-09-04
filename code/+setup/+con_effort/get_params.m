function outparams = get_params(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,...
    %
    % Any parameters shown here override the defaults set in
    % Classes/Params.m
        
    ii = 1;
    params(ii).name = sprintf('trial %i',ii); 
    params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
    params(ii).cmin = 0.001; % -1 for min(c(a,y)) of policy fn with no adj costs
    params(ii).cmax = 4; % -1 for max(c(a,y)) of policy fn with no adj costs
    params(ii).c_gcurv = 0.3;
    params(ii).bmax = 200;
    params(ii).bmin = - 0.0015 / (0.02/4);
    params(ii).b_gcurv_neg = 0.4;
    params(ii).rho = 0.007; % only relevant if runopts.IterateRho = 0 in master.m

    params(ii).rhoL = 0.01;

    params(ii).hdef = "cdot/c"; % "cdot" or "cdot/c"

    params(ii).chi0 = 0.5;
    params(ii).chi1 = 0.5;
    params(ii).chi2 = 1.5;
    params(ii).hbar = 0.1;
    
    params(ii).n_mpcsim = 2e4;
    params(ii).T_mpcsim = 1e3;
    
    params(ii).penalty1 = 1e3;
    params(ii).penalty2 = 2; % penalty for a < 0

    params(ii).nb = 200;
    params(ii).nb_neg = 20;
    params(ii).nb_KFE = 200;
    params(ii).nb_neg_KFE = 20;
    params(ii).b_gcurv_pos = 0.2;
    params(ii).nc = 200;
    params(ii).nc_KFE = 200;
    
    params(ii).Bequests = 1;
    
    % use h = cdot since cdot/c produces very small MPCs (~0.005)
    chi0s = [0.5,1,2]; % if < 1 then almost everyone responds to shock
    chi1s = [2,5,10]; % if < 5 then simulated mean MPCs are inaccurate
    chi2s = 1.5;
    
    ii = 2;
    for chi0 = chi0s
    for chi1 = chi1s
    for chi2 = chi2s
        params(ii).name = sprintf('trial %i',ii); 
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).cmin = 0.001;
        params(ii).cmax = 3;
        params(ii).c_gcurv = 0.25;
        params(ii).bmax = 150;
        params(ii).bmin = - 0.0015 / (0.02/4);
        params(ii).b_gcurv_neg = 0.7;
        params(ii).b_gcurv_pos = 0.08;
        params(ii).rho = 0.02; % only relevant if runopts.IterateRho = 0 in master.m
        params(ii).rhoL = 0.005;

        params(ii).hdef = "cdot"; % "cdot" or "cdot/c"
        
        params(ii).n_mpcsim = 5e4;
        params(ii).T_mpcsim = 1e3;

        params(ii).chi0 = chi0;
        params(ii).chi1 = chi1;
        params(ii).chi2 = chi2;
        params(ii).penalty1 = 1e3;
        params(ii).penalty2 = 2; % penalty for a < 0

        params(ii).nb = 200;
        params(ii).nb_neg = 40;
        params(ii).nb_KFE = 200;
        params(ii).nb_neg_KFE = 40;
        params(ii).nc = 300;
        params(ii).nc_KFE = 300;

        ii = ii + 1;
    end
    end
    end

    %% ------------------------------------------------------------------------
    % DO NOT CHANGE BELOW
    % -------------------------------------------------------------------------

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = setup.con_effort.ParamsConEffort(runopts,chosen_param);

end
