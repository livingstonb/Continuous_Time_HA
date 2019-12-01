classdef HACT_Model < handle
	properties
		p;

		grids_HJB;

		grids_KFE;

		income;

		hjb_solver;

		A_constructor;
	end

	methods
		function obj = HACT_Model(params, gridsHJB, gridsKFE, inc)
			obj.p = params;
			obj.grids_HJB = gridsHJB;
			obj.grids_KFE = gridsKFE;
			obj.income = inc;
		end

		function initialize(obj)
			import HACTLib.computation.HJBSolver
			import HACTLib.computation.HJBSolverSDU

			if obj.p.SDU
		        obj.hjb_solver = HJBSolverSDU(obj.p, obj.income, obj.p.hjb_options);
		    else
		        obj.hjb_solver = HJBSolver(obj.p, obj.income, obj.p.hjb_options);
			end

			returns_risk = obj.p.sigma_r > 0;
			obj.A_constructor_HJB = TransitionMatrixConstructor(obj.p, obj.income,...
				obj.grids_HJB, returns_risk);

			obj.kfe_solver = KFESolver(p, income, grdKFE, p.kfe_options);
		end

		function solve(obj)
			HJB = obj.solve_HJB();
			KFE = obj.solve_KFE();
		    
			%% --------------------------------------------------------------------
			% COMPUTE WEALTH
			% ---------------------------------------------------------------------
			wealth = compute_wealth(KFE.g, obj.grids_KFE);
			fprintf('    --- MEAN WEALTH = %f ---\n\n', wealth)
		end
	end

	methods (Access=protected)
		function HJB = solve_HJB(obj)
			import HACTLib.computation.make_initial_guess
			import HACTLib.computation.find_policies

			[Vguess, gguess] = make_initial_guess(p, obj.grids_HJB,...
					obj.grids_KFE, obj.income);

			fprintf('    --- Iterating over HJB ---\n')
		    dst = 1e5;
			for nn	= 1:obj.p.HJB_maxiters
				[HJB, V_deriv_risky_asset_nodrift] = find_policies(...
					p, income, grd, Vn);

			    % Construct transition matrix 
		        [A, stationary] = A_constructor_HJB.construct(HJB, Vn);

		        if obj.p.SDU
			    	risk_adj = {compute_risk_adjustment_for_nodrift_case(...
			    		p, grd, V_deriv_risky_asset_nodrift, stationary, Vn)};
			    else
			    	risk_adj = {};
			    end
		        
				% Update value function
				Vn1 = obj.hjb_solver.solve(A, HJB.u, Vn, risk_adj{:});

				% check for convergence
			    Vdiff = Vn1 - Vn;
			    Vn = Vn1;
			    dst = max(abs(Vdiff(:)));
		        if (nn==1) || (mod(nn,25)==0)
			    	fprintf('\tHJB iteration = %i, distance = %e\n', nn, dst);
		        end

		        if dst < obj.p.HJB_tol
			        fprintf('\tHJB converged after %i iterations\n', nn);
				 	break
				end
				check_if_not_converging(dst, nn);
		    end

		    if (nn >= p.HJB_maxiters)
		        error("HJB didn't converge");
		    end
		    
		    % Store value function and policies on both grids
		    HJB.Vn = Vn;
		end

		%% --------------------------------------------------------------------
	    % SOLVE KFE
		% ---------------------------------------------------------------------
		function solve_KFE(obj, HJB)
			KFE.Vn = reshape(interp_decision * HJB.Vn(:),...
		    	[obj.p.nb_KFE, obj.p.na_KFE, obj.p.nz, obj.income.ny]);
		    KFE = find_policies(obj.p, obj.income, obj.grids_KFE, KFE.Vn);

			import HACTLib.computations.TransitionMatrixConstructor

			% True if returns should be treated as risky in the KFE
			returns_risk = (p.sigma_r > 0) && (p.retrisk_KFE == 1);
		    A_constructor_kfe = TransitionMatrixConstructor(p, income, grdKFE, returns_risk);
		    Au = A_constructor_kfe.construct(KFE, KFE.Vn);
		    
		    
		    KFE.g = kfe_solver.solve(Au);
		end
	end
end

function risk_adj = compute_risk_adjustment_for_nodrift_case(...
	p, grd, V_deriv_risky_asset_nodrift, stationary, Vn)
	% computes the adjustment term when returns are risky
	% and there is no drift in the risky assets for some
	% asset > 0 cases
	
	% Parameters
	% ---------
	% p : a Params object
	%
	% grd : a Grid object
	%
	% V_deriv_risky_asset_nodrift : approximation of the
    %	first derivative of V for the case with no drift
    %	in the risky asset
    %
    % stationary : boolean mask to indicate states where
    %	there is no drift in risky asset
    %
    % Vn : the value function over states
    %
    % Returns
    % -------
    % risk_adj : an array of shape (nb, na, nz, ny) with
    % 	zeros everywhere except for states where risky
    %	asset drift is zero, where the array contains
    %	the term with Va^2 

	if ~isempty(stationary)
		% there are states with neither backward nor forward drift,
		% need to compute additional term for risk
		if p.invies == 1
			risk_adj = (1-p.riskaver) * V_deriv_risky_asset_nodrift .^ 2;
		else
			risk_adj = V_deriv_risky_asset_nodrift .^ 2 ./ Vn * (p.invies - p.riskaver) / (1-p.invies);
		end

		risk_adj = risk_adj .* (grd.a.matrix * p.sigma_r) .^ 2 / 2;
		risk_adj(~stationary) = 0;
	else
		risk_adj = [];
	end
end

function check_if_not_converging(dst, nn)
	if dst>10 && nn>500
	 	% Not going to converge, throw exception
	 	msgID = 'HACTLib:HACT_Model:NotConverging';
	    msg = 'The HJB does not appear to be converging';
	    HJBException = MException(msgID,msg);
	    throw(HJBException)
	end
end