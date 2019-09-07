classdef Params < handle
    % This class stores the parameters of the model and contains methods to
    % adjust them for frequency and other factors. See properties for the
    % default values.
    
    properties (SetAccess=protected)
        % Run options
        IterateRho;
        ComputeMPCS;
        SimulateMPCS;
        ComputeMPCS_news;
        SimulateMPCS_news;
        Bequests = 0;
        ResetIncomeUponDeath = 0; % WARNING: keep = 0, may not be up to date
        
        % Identifier
        name = 'unnamed';
        param_index;

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

        % income grid size
        ny;

        % other heterogeneity
        nz;

        min_grid_spacing = 0.001;
        

	    % ------------ returns -----------------------------
		r_b = 0.02 / 4; %0.005; % steady-state savings rate, liquid assets
		perfectannuities = 0;

		% ------------ preferences -------------------------
		riskaver = 1; % risk aversion
		deathrate = 1 /200; % death rate (quarterly)
		rho = 0.015; % discount factor (quarterly)
		rhoL; % starting point to find rho lower bounds
		rhos;
		rho_grid;

		% ------------ income process ---------------------
		DirIncomeProcess = 'input/IncomeGrids/continuous_a';

		% ------------ taxes ------------------------------
		transfer = 0; % transfer to households 
		wagetax = 0; % tax rate on wage income

		% ------------ targets ----------------------------
		targetAY = 3.5; % target for mean assets / mean annual income

		% ------------ approximation parameters -----------
		% HJB loop
		maxit_HJB 		= 600; % maximal allowable number of HJB iterations
		crit_HJB 		= 1e-6; % critical value
		delta_HJB 		= 1e6; % step size (large steps work even though method is partially explicit)

		% Howard improvement step in HJB loop
		maxit_HIS 		= 10; % total number of Howard improvement steps
		start_HIS 		= 2; % when in HJB loop do Howard improvement steps begin?
		crit_HIS 		= 1e-5; % critical value

		% KFE loop
		maxit_KFE 		= 1e4; % maximal allowable number of KFE iterations
		crit_KFE 		= 1e-7; % critical value
		delta_KFE 		= 1e6; %1e6; % step size

		% Outer assets-income ratio grid
		maxit_AY 		= 100; % maximal allowable number of loops over capital-labor ratio
		crit_AY			= 1e-7; % critical value

        % Step-size for Feynman-Kac formula
        delta_mpc = 0.025;
        
        % ----------- statistics variables -------------------------------
        epsilon_HtM = [0 0.005 0.01 0.02 0.05 0.1 0.15]; % for looking at fraction HtM
        wpercentiles = [10 25 50 90 99 99.9];
        mpc_shocks = [-1e-5 -0.01 -0.1 1e-5 0.01 0.1];
        
        T_mpcsim = 1e3;
        n_mpcsim = 1e5;
	end


    methods
        function obj = Params(runopts,params)
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
            
            obj.param_index = runopts.param_index;
            obj.IterateRho = runopts.IterateRho;
            obj.ComputeMPCS = runopts.ComputeMPCS;
            obj.SimulateMPCS = runopts.SimulateMPCS;
            obj.ComputeMPCS_news = runopts.ComputeMPCS_news;
            obj.SimulateMPCS_news = runopts.SimulateMPCS_news;
            obj.tempdirec = runopts.temp;

            obj.rhos = obj.rho + obj.rho_grid;
        end

        function update_ny(obj,ny)
        	obj.ny = ny;
        end
    end
end
