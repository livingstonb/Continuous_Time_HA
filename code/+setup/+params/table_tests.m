function outparams = table_tests(param_opts, param_index)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,.

    import HACTLib.aux.set_shared_fields

    shared_params = param_opts;
    shared_params.chi0 = 0;
    shared_params.a_lb = 0.25;
    shared_params.OneAsset = 0;
    shared_params.calibration_vars = {'rho'};
    shared_params.calibration_stats = {'totw'};
    shared_params.calibration_targets = [3.5];
    shared_params.calibration_bounds = {[0.001, 0.05]};
    shared_params.income_dir = 'continuous_b';
    shared_params.nb_KFE = 120;
    shared_params.nb = 120;
    shared_params.na = 120;
    shared_params.na_KFE = 120;
    shared_params.b_gcurv_pos = 0.2;
    shared_params.a_gcurv = 0.2;
    shared_params.bmax = 100; %5;
    shared_params.amax = 200; %20;
%     shared_params.min_grid_spacing = 1e-5;

    %% adjusing (A/Y,B/Y), fixed adj costs
    raBASELINE = 0.0190643216;

    ii = 1;
    if param_index < 100
	    % ras = [raBASELINE, 0.0129226573755814, 0.0304148974027774];
	    ras = raBASELINE;
	    % calibrate r_a to match B/Y in each case
	    for ra = ras
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('contB, fixed_adj_costs, r_a=%f',ra);
	        params{ii}.chi1 = 0.15;
	        params{ii}.chi2 = 0.25;
	        params{ii}.rho = 0.015440584992491;
	        params{ii}.r_a = ra;
	        params{ii}.calibration_vars = {'rho', 'r_a'};
	        params{ii}.calibration_stats = {'totw', 'liqw'};
	        params{ii}.calibration_targets = [3.5, 0.5];
	        params{ii}.calibration_bounds = {[0,001, 0.05], [0.006, 0.05]};
	        
	        ii = ii + 1;
	    end
	elseif param_index < 200  
	    %% varying chi1 and chi2
	    for chi1 = [0.1,0.4]
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('chi1=%f',chi1);
	        params{ii}.chi1 = chi1;
	        params{ii}.chi2 = 0.25;
	        params{ii}.r_a = raBASELINE;
	        
	        ii = ii + 1;
	    end
	    
	    for chi2 = [0.1,0.15,0.5,1]
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('chi2=%f',chi2);
	        params{ii}.chi1 = 0.15;
	        params{ii}.chi2 = chi2;
	        params{ii}.r_a = raBASELINE;
	        
	        ii = ii + 1;
	    end
	    
	    %% varying r_a and r_b
	    
	    ras_annual = 0.05:0.01:0.1;
	    for i = 1:6
	        ras_quarterly(i) = (1+ras_annual(i)).^(1/4) - 1;
	    end
	elseif param_index < 300
		   
	    for ra = ras_quarterly
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('annual r_a=%f', ras_annual(ii));
	        params{ii}.chi1 = 0.15;
	        params{ii}.chi2 = 0.25;
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
	        params{ii}.chi1 = 0.15;
	        params{ii}.chi2 = 0.25;
	        params{ii}.r_b = rb;
	        params{ii}.r_a = 0.0190643216;
	        
	        ii = ii + 1;
	    end
	elseif param_index < 400
	    %% recalibrate (rho, r_a) to match B/Y=0.5 with cont_a
	    ras = 0.0173096090154433;
	    for ra = ras
	    	params{ii} = shared_params;
	        params{ii}.name = 'cont_a, recalibrated';
	        params{ii}.income_dir = 'continuous_a';
	        params{ii}.chi1 = 0.15;
	        params{ii}.chi2 = 0.25;
	        params{ii}.r_a = ra;

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
	    params{ii}.chi1 = 0.15;
	    params{ii}.chi2 = 0.25;
	    params{ii}.r_a = 0.01906432159025;
	    params{ii}.rho = 0.015440584980253;
	elseif param_index < 600
	    %% risk aversion coeff tests
	    for riskaver = [0.5,2,4,6]
	    	params{ii} = shared_params;
	        params{ii}.name = sprintf('risk aversion = %f',riskaver); 
	        params{ii}.chi1 = 0.15;
	        params{ii}.chi2 = 0.25;
	        params{ii}.rho = 0.005;
	        params{ii}.r_a = raBASELINE;
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
	        params{ii}.name = sprintf('rho het, width = %f',w);
	        params{ii}.chi1 = 0.15;
	        params{ii}.chi2 = 0.25;
	        params{ii}.rho = 0.005;
	        params{ii}.r_a = raBASELINE;
	        params{ii}.rho_grid = [-w, 0, w];
	        ii = ii + 1;
	    end
	end

    %% DO NOT CHANGE THIS SECTION
    % Use param_opts.param_index to choose which specification to select
    chosen_param = params{mod(param_index, 100)};
end
