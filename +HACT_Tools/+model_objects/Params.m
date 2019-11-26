classdef Params < handle
    % This class stores the parameters of the model and contains methods to
    % adjust them for frequency and other factors. See properties for the
    % default values.
    
    properties (SetAccess=protected)
        % Run options
        ComputeMPCS;
        SimulateMPCS;
        ComputeMPCS_news;
        SimulateMPCS_news;
        Bequests = 0;
        ResetIncomeUponDeath = 0; % WARNING: keep = 0, may not be up to date
        SaveResults = 1;
        OneAsset = 1;
        DealWithSpecialCase;
        NoRisk = 1;
        
        % Identifier
        name = 'unnamed';
        param_index;

        maindir;
        income_dir;
        tempdirec;

        % ------------ grid parameters ---------------------
        % liquid grid parameters
		bmin = 0;
		b_soft_constraint = 0;
		bmax = 50;
		b_gcurv_pos = 0.2;
		b_gcurv_neg = 0.4;
        nb = 150;
		nb_pos; % default is set to nb
        nb_neg;
		nb_KFE = 150;
		nb_pos_KFE; % default is set to nb_KFE
        nb_neg_KFE;

        % illiquid grid parameters
        na = 120;
        amin = 0;
        amax = 200;
        a_gcurv = 0.2;
        na_KFE = 120;

        % income grid size
        ny;

        % other heterogeneity
        nz;

        min_grid_spacing = 0.001;
        

	    % ------------ returns -----------------------------
		r_b = 0.02 / 4; %0.005; % steady-state savings rate, liquid assets
        r_a = 0.06/4; % illiquid return
		perfectannuities = 0;

        % borrowing
        borrwedge = 0.08/4; % steady-state borrowing wedge, liquid assets
        r_b_borr; % steady-state borrowing rate, liquid assets

        % Rate of return risk
        sigma_r = 0; % turn off by default
        retrisk_KFE = 0; % by default, don't include risk in KFE when sigma_r > 0

		% ------------ preferences -------------------------
        SDU = 0;
		riskaver = 1; % risk aversion, can be a row vector
        % Stochastic differential utility
        invies = 1;
		riskaver_fulldim;
		deathrate = 1 /200; % death rate (quarterly)
		rho = 0.015; % discount factor (quarterly)
		rhoL; % starting point to find rho lower bounds
		rhos;
		rho_grid;

		% ------------ taxes ------------------------------
		transfer = 0; % transfer to households 
		wagetax = 0; % tax rate on wage income

		% ------------ targets ----------------------------
		targetAY = 3.5; % target for mean assets / mean annual income

		% ------------ approximation parameters -----------
		% HJB loop
		hjb_options;
		HJB_maxiters = 2000; % maximal allowable number of HJB iterations
		HJB_tol = 1e-8; % critical value
		HJB_delta = 1e6; % step size
        HJB_implicit = false;

		% Howard improvement step in HJB loop
		HIS_maxiters = 10; % total number of Howard improvement steps
		HIS_start = 2; % when in HJB loop do Howard improvement steps begin?
		HIS_tol= 1e-5; % critical value

		% KFE loop
		kfe_options;
		KFE_maxiters = 1e4; % maximal allowable number of KFE iterations
		KFE_tol = 1e-8; % critical value
		KFE_delta = 1e6; %1e6; % step size
        KFE_iterative = true;

		% Outer assets-income ratio grid
		maxit_AY 		= 100; % maximal allowable number of loops over capital-labor ratio
		crit_AY			= 1e-7; % critical value

        % Step-size for Feynman-Kac formula
        delta_mpc = 0.025;

        % marginal product of labor/capital
        endogenousLabor = 0;
        labor_disutility = 1;
        frisch = 0.5;
        MPL = 1;
        MPK = 1;
        
        % ----------- statistics variables -------------------------------
        epsilon_HtM = [0 0.005 0.01 0.02 0.05 0.1 0.15]; % for looking at fraction HtM
        wpercentiles = [10 25 50 90 99 99.9];
        mpc_shocks = [-1e-5 -0.01 -0.1 1e-5 0.01 0.1];
        decomp_thresholds = [0 0.01 0.05];
        
        T_mpcsim = 500;
        n_mpcsim = 1e5;

        % ----------- adjustment costs -----------------------------------
        chi0 = 0.070046;
        chi1 = 0.535791;
        chi2 = 0.971050;
        a_lb = 0.2407;

        % deposit share     
        directdeposit = 0;
	end

    methods
        function obj = Params(runopts,params)

        	import HACT_Tools.options.HJBOptions
        	import HACT_Tools.options.KFEOptions

			 % allow some properties to be set by params structure
            if nargin > 1
                if isstruct(params)
                	pfields = fields(params);
		            if ~all(ismember(pfields,properties(obj)))
		            	indices = find(~ismember(pfields,properties(obj)));
			    		msg = sprintf('Some field(s) in the input parameters is invalid:\n');
			    		for ii = indices
			    			msg = [msg sprintf('\t%s',pfields{ii})];
			    		end
			    		error(msg)
                    elseif numel(params) > 0
		                for ip = 1:numel(pfields)
                            if ~isempty(params.(pfields{ip}))
                                obj.(pfields{ip}) = params.(pfields{ip});
                            end
		                end
		            end
			    else
		         	error('Invalid argument type passed to constructor')
                end
            end

            % check for other heterogeneity
            if (numel(obj.rho_grid)>1) && (numel(obj.riskaver)>1)
            	error('Cannot have both rho and riskaver heterogeneity')
            else
                obj.nz = max(numel(obj.rho_grid),numel(obj.riskaver));
            end

             obj.kfe_options = KFEOptions(...
				'maxiters', obj.KFE_maxiters,...
				'tol', obj.KFE_tol,...
				'delta', obj.KFE_delta,...
				'iterative', obj.KFE_iterative);

            obj.hjb_options = HJBOptions(...
				'delta', obj.HJB_delta,...
				'implicit', obj.HJB_implicit,...
				'HIS_maxiters', obj.HIS_maxiters,...
				'HIS_tol', obj.HIS_tol,...
				'HIS_start', obj.HIS_start);
            
            obj.param_index = runopts.param_index;
            obj.ComputeMPCS = runopts.ComputeMPCS;
            obj.SimulateMPCS = runopts.SimulateMPCS;
            obj.ComputeMPCS_news = runopts.ComputeMPCS_news;
            obj.SimulateMPCS_news = runopts.SimulateMPCS_news;
            obj.maindir = runopts.direc;
            obj.tempdirec = runopts.temp;

            obj.rhos = obj.rho + obj.rho_grid;
            obj.riskaver_fulldim = reshape(obj.riskaver,[1 1 numel(obj.riskaver) 1]);

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
        
        function obj = set(obj, field, new_val, quiet)
        	% Sets the value of a parameter.
        	%
        	% Inputs
        	% ------
        	%
        	% field : A string containing the parameter
        	%	name.
        	%
        	% new_val : The desired value of the parameter.
        	%
        	% quiet : An optional argument that, when it
        	%	evaluates to true, suppresses printing
        	%	to the screen.

        	field = char(field);

        	KFE_option_passed = false;
        	if numel(field) > 4
        		if strcmp(field(1:3), 'KFE')
        			KFE_option_passed = true;
        		end
        	end

        	if strcmp(field, 'rho')
        		obj.rhos = new_val + obj.rho_grid;
        	elseif KFE_option_passed
        		assert(isprop(obj.kfe_options, field(5:end)), "Invalid KFE option");
        		obj.kfe_options.set(field(5:end), new_val);
    		elseif ~isprop(obj, field)
    			error("Requested field is not an attribute of Params.");
        	end

            obj.(field) = new_val;

            if ~exist('quiet', 'var')
            	quiet = false;
            end

            if ~quiet
            	disp(strcat(field, sprintf(" has been reset to %.9f", new_val)));
            end
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
