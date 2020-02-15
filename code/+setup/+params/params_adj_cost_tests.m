function outparams = params_adj_cost_tests(runopts)

	chi1s = [0.01, 0.05, 0.1, 0.15, 0.2];
    chi2s = [0.1, 25, 0.5];

    %% --------------------------------------------------------------------
    % CREATE PARAMETERIZATIONS HERE
    % ---------------------------------------------------------------------
    ii = 1;
    for chi1 = chi1s
    	for chi2 = chi2s
		    params(ii).name = sprintf('chi1=%g, chi2=%g',chi1, chi2); 
		    params(ii).OneAsset = 0;
		    params(ii).income_dir = 'continuous_b';
		    params(ii).chi0 = 0;
		    params(ii).chi1 = chi1;
		    params(ii).chi2 = chi2;
		    params(ii).a_lb = 0.25;
		    params(ii).rho = 0.015440584992491;
		    params(ii).r_a = 0.02;
		    ii = ii + 1;
		end
	end

    %% DO NOT CHANGE BELOW

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACTLib.model_objects.Params(runopts,chosen_param);
end