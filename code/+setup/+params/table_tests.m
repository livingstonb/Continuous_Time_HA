function outparams = table_tests(param_opts, param_index)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,.

    import HACTLib.aux.set_shared_fields

    shocks = [-1, -500, -5000, 1, 500, 5000];
    dollars = 72000;

    shared_params = param_opts;
    shared_params.mpc_shocks = shocks / dollars;
    shared_params.numeraire_in_dollars = dollars;
    shared_params.chi0 = 0;
    shared_params.a_lb = 0.25;
    shared_params.chi1 = 0.15;
	shared_params.chi2 = 1;
    shared_params.OneAsset = 0;
    shared_params.calibration_vars = {'rho'};
    shared_params.calibration_stats = {'totw'};
    shared_params.calibration_targets = [3.5];
    shared_params.calibration_bounds = {[0.001, 0.05]};
    shared_params.income_dir = 'continuous_b';
    shared_params.nb_KFE = 80;
    shared_params.nb = 80;
    shared_params.na = 80;
    shared_params.na_KFE = 80;
    shared_params.b_gcurv_pos = 0.2;
    shared_params.a_gcurv = 0.2;
    shared_params.a_glinear = 0.02;
    shared_params.b_glinear = 0.01;
    shared_params.bmax = 20;
    shared_params.amax = 200;
    shared_params.rho = 0.0144188;
    shared_params.r_a = 0.0188616;

    ii = 1;
    if param_index < 100
	    % ras = [raBASELINE, 0.0129226573755814, 0.0304148974027774];
	    % calibrate r_a to match B/Y in each case
    	params{ii} = shared_params;
        params{ii}.name = 'baseline';
        params{ii}.r_a = 0.015;
        params{ii}.rho = 0.012;
        params{ii}.calibration_vars = {'rho', 'r_a'};
        params{ii}.calibration_stats = {'totw', 'liqw'};
        params{ii}.calibration_targets = [3.5, 0.5];
        params{ii}.calibration_bounds = {[0,001, 0.05], [0.006, 0.05]};
        
        ii = ii + 1;
	elseif param_index < 200  
	    %% varying chi1 and chi2
	    for chi1 = [0.1,0.4]
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('chi1=%f',chi1);
	        params{ii}.chi1 = chi1;
	        
	        ii = ii + 1;
	    end
	    
	    for chi2 = [0.1,0.15,0.5,1]
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('chi2=%f',chi2);
	        params{ii}.chi2 = chi2;
	        
	        ii = ii + 1;
	    end
	    
	elseif param_index < 300

		%% varying r_a and r_b
	    
	    ras_annual = 0.05:0.01:0.1;
	    for i = 1:6
	        ras_quarterly(i) = (1+ras_annual(i)).^(1/4) - 1;
	    end
		   
	    for ra = ras_quarterly
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('annual r_a=%f', ras_annual(ii));
	        params{ii}.r_a = ra;
	        
	        ii = ii + 1;
	    end
	    
	    rbs_annual = [0,0.05];
	    for i = 1:2
	        rbs_quarterly(i) = (1+rbs_annual(i)).^(1/4) - 1;
	    end
	    
	    for rb = rbs_quarterly
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('annual r_b=%f',(1+rb)^4-1);
	        params{ii}.r_b = rb;
	        
	        ii = ii + 1;
	    end
	elseif param_index < 400
	    for ra = ras
	    	params{ii} = shared_params;
	        params{ii}.name = 'cont_a, recalibrated';
	        params{ii}.income_dir = 'continuous_a';
	        params{ii}.rho = 0.0127986;
	        params{ii}.r_a = 0.0176018;

	        params{ii}.calibration_vars = {'rho', 'r_a'};
	        params{ii}.calibration_stats = {'totw', 'liqw'};
	        params{ii}.calibration_targets = [3.5, 0.5];
	        params{ii}.calibration_bounds = {[0,001, 0.05], [0.006, 0.05]};
	        
	        ii = ii + 1;
	    end
	elseif param_index < 500
	    
	    %% run cont_a without re-calibrating (rho,r_a)
	    params{ii} = shared_params;
	    params{ii}.name = sprintf('conta, no recalibration'); 
	    params{ii}.income_dir = 'continuous_a';

	elseif param_index < 600
	    %% risk aversion coeff tests
	    for riskaver = [0.5,2,4,6]
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('risk aversion = %f',riskaver); 
	        params{ii}.riskaver = riskaver;
	        
	        if riskaver == 4
	            params{ii}.KFE_delta = 0.2;
	            params{ii}.KFE_maxiters = 25000;
	        elseif riskaver == 6
	            params{ii}.KFE_delta = 0.01;
	        end
	        ii = ii + 1;
	    end
	elseif param_index < 700

	    %% rho heterogeneity
	    for w = [0.001, 0.005]
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('rho het, spacing = %f',w);
	        params{ii}.rho = 0.005;
	        params{ii}.rho_grid = [-w, 0, w];
	        ii = ii + 1;
	    end
	end

    %% DO NOT CHANGE THIS SECTION
    % Use param_opts.param_index to choose which specification to select
    outparams = params{mod(param_index, 100)};
end
