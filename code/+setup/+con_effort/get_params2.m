function outparams = get_params2(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,...
    %
    % Any parameters shown here override the defaults set in
    % Classes/Params.m
    
    targetAYs = [0.5,1,2];
    chi0s = [1,2,5];
    hbars = [5,10];
    hdefs = {'cdot','cdot/c'};
    
    ii = 1;
    for targetAY = targetAYs
    for chi0 = chi0s
    for hbar = hbars
    for hdef = hdefs
        params(ii).name = sprintf('wealth=%2.1f, chi0=%i, hbar=%i, hdef=%s',targetAY,chi0,hbar,hdef{1}); 
        params(ii).DirIncomeProcess = 'input/IncomeGrids/continuous_b';
        params(ii).cmin = 0.001;
        params(ii).cmax = 3;
        params(ii).c_gcurv = 0.25;
        params(ii).bmax = 150;
        params(ii).bmin = - 0.0015 / (0.02/4);
        params(ii).b_gcurv_neg = 0.7;
        params(ii).b_gcurv_pos = 0.08;
        params(ii).rho = 0.02; % only relevant if runopts.IterateRho = 0 in master.m
        params(ii).rhoL = 0.03;

        params(ii).hdef = hdef{1}; % "cdot" or "cdot/c"

        params(ii).targetAY = targetAY;
        
        params(ii).n_mpcsim = 5e4;
        params(ii).T_mpcsim = 1e3;

        params(ii).chi0 = chi0;
        params(ii).chi1 = 0;
        params(ii).hbar = hbar;
        params(ii).penalty1 = 1e3;
        params(ii).penalty2 = 2; % penalty for a < 0

        params(ii).nb = 200;
        params(ii).nb_pos = 170;
        params(ii).nb_KFE = 200;
        params(ii).nb_pos_KFE = 170;
        params(ii).nc = 300;
        params(ii).nc_KFE = 300;

        params(ii).Bequests = 1;

        ii = ii + 1;
    end
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
