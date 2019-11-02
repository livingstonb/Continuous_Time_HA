classdef ParamsTwoAsset < setup.Params

	properties (SetAccess={?setup.Params})
		ModelType = 'NoEffortCost';
		OneAsset = 1;
		DealWithSpecialCase;
        NoRisk = 1;

        % illiquid grid parameters
        na = 120;
        amin = 0;
        amax = 200;
        a_gcurv = 0.2;
        na_KFE = 120;

        % illiquid return
        r_a = 0.06/4;

        % borrowing
        borrwedge = 0.08/4; % steady-state borrowing wedge, liquid assets
		r_b_borr; % steady-state borrowing rate, liquid assets

        % adjustment cost function
		chi0 = 0.070046;
		chi1 = 0.535791;
		chi2 = 0.971050;
		a_lb = 0.2407;

		% Deposit share		
		directdeposit = 0;

		decomp_thresholds = [0 0.01 0.05];

        % Stochastic differential utility
        invies = 1;

        % Rate of return risk
        sigma_r = 0; % turn off by default
        retrisk_KFE = 1; % by default, include risk in KFE when sigma_r > 0
	end

	methods
		function obj = ParamsTwoAsset(runopts,params)
			obj = obj@setup.Params(runopts,params);

			obj.SimulateMPCS = runopts.SimulateMPCS;
			obj.DealWithSpecialCase = runopts.DealWithSpecialCase;

			% adjust interest rates
            obj.r_b_borr = obj.r_b + obj.borrwedge;

            % adjust if oneasset is selected
            if obj.OneAsset == 1
                obj.na = 2;
                obj.na_KFE = 2;
                obj.nb = 500;
                obj.nb_KFE = 400;
                obj.b_gcurv_pos = 0.2;
                obj.b_gcurv_neg = 0.2;
                obj.a_gcurv = 0.2;
                obj.bmax = 100;
                obj.chi0 = 1e8;
                obj.r_a = obj.r_b;
            end

            if runopts.fast == 1
                obj.nb = 16;
                obj.nb_pos = 16;
                obj.nb_neg = 0;
                obj.nb_KFE = 15;
                obj.nb_pos_KFE = 15;
                obj.nb_neg_KFE = 0;

                % illiquid grid parameters
                if obj.OneAsset == 0
                    obj.na = 14;
                    obj.na_KFE = 13;
                end
                
            	obj.n_mpcsim = 100;
            	obj.T_mpcsim = 1e3;
            end
            
            % Set default grid sizes
            if isempty(obj.nb_pos) == 1
                obj.nb_pos = obj.nb;
            end
            if isempty(obj.nb_pos_KFE) == 1
                obj.nb_pos_KFE = obj.nb_KFE;
            end
            
            obj.nb_neg = obj.nb - obj.nb_pos;
            obj.nb_neg_KFE = obj.nb_KFE - obj.nb_pos_KFE;
            
            if obj.OneAsset == 1
                obj.rhoL = obj.r_a + obj.perfectannuities*obj.deathrate - obj.deathrate;
                obj.rhoL = obj.rhoL + 2.5e-3;
            elseif ~ismember('rhoL',fields(params))
                % two-asset case but rhoL was not set in parameters file
                obj.rhoL = 0.023;
            end
		end

        function reset_returns(obj, r_b, r_a)
            obj.r_b = r_b;
            obj.r_a = r_a;
        end
		
		function print(obj)
        	fprintf('\n\nSelected parameterization %i:\n',num2str(obj.param_index)) 
            fprintf('%s\n\n',obj.name)

            fprintf('Chosen parameters were...\n\n')

            fprintf('\tb_soft_constraint = %f\n',obj.b_soft_constraint)
            fprintf('\tbmin = %f\n',obj.bmin)
			fprintf('\tbmax = %f\n',obj.bmax)
			fprintf('\tb_gcurv_pos = %f\n',obj.b_gcurv_pos)
			fprintf('\tb_gcurv_neg = %f\n',obj.b_gcurv_neg)
			fprintf('\tnb = %i\n',obj.nb)
			fprintf('\tnb_pos = %i\n',obj.nb_pos)
			fprintf('\tnb_KFE = %i\n',obj.nb_KFE)
			fprintf('\tnb_pos_KFE = %i\n',obj.nb_pos_KFE)
			fprintf('\n')

			if obj.OneAsset == 0
				fprintf('\tna = %i\n',obj.na)
				fprintf('\tna_KFE = %i\n',obj.na_KFE)
				fprintf('\tamin = %f\n',obj.amin)
				fprintf('\tamax = %f\n',obj.amax)
				fprintf('\ta_gcurv = %f\n',obj.a_gcurv)
				fprintf('\n')

				fprintf('\tchi0 = %f\n',obj.chi0)
				fprintf('\tchi1 = %f\n',obj.chi1)
				fprintf('\tchi2 = %f\n',obj.chi2)
				fprintf('\ta_lb = %f\n',obj.a_lb)
				fprintf('\n')
			end

			fprintf('\tBequests = %i\n',obj.Bequests)
			fprintf('\tResetIncomeUponDeath %i\n',obj.ResetIncomeUponDeath)
			fprintf('\n')

			fprintf('\tr_b = %f\n',obj.r_b)
			fprintf('\tr_b_borr = %f\n',obj.r_b_borr)
			fprintf('\tr_a = %f\n',obj.r_a)
			fprintf('\tperfectannuities = %f\n',obj.perfectannuities)
			fprintf('\triskaver = %f\n',obj.riskaver)
			fprintf('\tdeathrate = %f\n',obj.deathrate)
            fprintf('\n\n')
        end
	end
end