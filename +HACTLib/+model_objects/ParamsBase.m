classdef ParamsBase < handle
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
        
        % Name for the parameterization
        name = 'unnamed';

        % Numeric identifier for the parameterization
        param_index;

        % Path to the root directory
        main_dir;

        % Path to the income variables
        income_dir;

        % Path to the temp directory
        temp_dir;

        %% -------------------------------------------
		% Liquid Asset Grid Parameters
		% --------------------------------------------
        % Min value of liquid assets
		bmin = 0;

		% Minimum value of liquid assets before
		% borrowing wedge kicks in
		b_soft_constraint = 0;

		% Max value for liquid assets
		bmax = 50;

		% Curvature of positive section of liquid asset,
		% lower value implies more curvature
		b_gcurv_pos = 0.2;

		% Curvature of negative section of liquid asset,
		% lower value implies more curvature
		b_gcurv_neg = 0.4;

		% Number of points on liquid asset grid
        nb = 150;
        
        % Number of points on positive section of liquid
        % asset grid, defaults to nb
		nb_pos;

		% Number of points on negative section of liquid
		% asset grid, defaults to 0
        nb_neg;

        % Number of points on liquid asset grid for the KFE
		nb_KFE = 150;

		% Number of points on positive section of liquid
		% asset grid for the KFE, defaults to nb_KFE
		nb_pos_KFE;

		% Number of points on negative section of liquid
		% asset grid for the KFE, defaults to 0
        nb_neg_KFE;

        %% -------------------------------------------
		% Illiquid Asset Grid Parameters
		% --------------------------------------------
		% Min value of illiquid assets
		amin = 0;

		% Max value of illiquid assets
        amax = 200;

        % Number of points on the illiquid asset grid
        na = 120;
        
        % Curvature of positive section of illiquid asset,
		% lower value implies more curvature
        a_gcurv = 0.2;

        % Number of points on illiquid asset grid for
        % the KFE
        na_KFE = 120;

        %% -------------------------------------------
		% Other Parameters
		% --------------------------------------------
        % Number of income states
        ny;

        % Number of states in z-dimension
        nz;

        % Minimum spacing of the asset grids close to the
        % constraints
        min_grid_spacing = 0.001;
        
	    %% -------------------------------------------
		% Returns
		% --------------------------------------------
		% Liquid returns, at a quarterly rate
		r_b = 0.02 / 4;

		% Illiquid returns, at a quarterly rate
        r_a = 0.06/4;

        % Annuities, 0 turns off annuities
		perfectannuities = 0;

        % Wedge between borrowing rate and r_b on the
        % liquid asset, only relevant with bmin < 0
        borrwedge = 0.08/4;

        % Liquid asset borrowing rate, set according
        % to borrwedge
        r_b_borr;

        % Standard deviation of returns risk
        sigma_r = 0;

        % Boolean indicator, set to true to include
        % returns risk in the KFE
        retrisk_KFE = false;

		%% -------------------------------------------
		% Preferences
		% --------------------------------------------

		% Stochastic differential utility
        SDU = false;

        % Coefficient of risk aversion, can be set
        % to a row vector to accomodate RA heterogeneity
		riskaver = 1;

        % Inverse of the IES, only relevant if using
        % SDU preferences
        invies = 1;

        % Array of risk aversion coefficients, set
        % automatically
		riskaver_fulldim;

		% Transition into death, at a quarterly rate
		deathrate = 1 /200;

		% Time discount factor, at a quarterly rate
		rho = 0.015;

		% Grid for rho values if there is rho
		% heterogeneity, automatically set to
		% rho + rho_grid
		rhos;

		% Grid to accomodate rho heterogeneity
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
end