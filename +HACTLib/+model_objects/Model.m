classdef Model < handle
	properties
		p;

		grids_HJB;

		grids_KFE;

		income;

		hjb_solver;

		kfe_solver;

		A_constructor_HJB;
	end

	methods
		function obj = Model(params, gridsHJB, gridsKFE, inc)
			obj.p = params;
			obj.grids_HJB = gridsHJB;
			obj.grids_KFE = gridsKFE;
			obj.income = inc;
		end

		function initialize(obj)
			import HACTLib.computation.TransitionMatrixConstructor
			import HACTLib.computation.HJBSolver
			import HACTLib.computation.KFESolver

			obj.hjb_solver = HJBSolver(obj.p, obj.income, obj.p.hjb_options);

			returns_risk = obj.p.sigma_r > 0;
			obj.A_constructor_HJB = TransitionMatrixConstructor(obj.p, obj.income,...
				obj.grids_HJB, returns_risk);

			obj.kfe_solver = KFESolver(obj.p, obj.income,...
				obj.grids_KFE, obj.p.kfe_options);
		end

		function [HJB, KFE, Au] = solve(obj)
			if obj.p.endogenous_labor
				% Find hours policy function at bmin over wage values.
				import HACTLib.model_objects.CRRA
				import HACTLib.model_objects.Frisch

				hours_bc = 0.5 * ones(obj.income.ny, 1);
				% hours_bc_last = hours_bc;
				for ih = 1:obj.p.HOURS_maxiters
					con = (obj.p.r_b+obj.p.deathrate*obj.p.perfectannuities)...
						* obj.grids_HJB.b.vec(1) + (1-obj.p.directdeposit-obj.p.wagetax)...
						* hours_bc .* obj.income.y.vec + obj.p.transfer;
					u1 = CRRA.marginal_utility(con, obj.p.invies);
					v1 = (1-obj.p.directdeposit-obj.p.wagetax)...
						* hours_bc .* obj.income.y.vec .* u1;
					hours_bc = Frisch.inv_marginal_disutility(v1,...
						obj.p.labor_disutility, obj.p.frisch);
					% hours_bc = min(hours_bc, 1);
					% max(abs(hours_bc(:)-hours_bc_last(:)))
					% hours_bc_last = hours_bc;
				end
				hours_bc = reshape(hours_bc, [1, 1, 1, obj.income.ny]);
				hours_bc_HJB = repmat(hours_bc, [1, obj.p.na, obj.p.nz, 1]);
				hours_bc_KFE = repmat(hours_bc, [1, obj.p.na_KFE, obj.p.nz, 1]);
			else
				hours_bc_HJB = [];
				hours_bc_KFE = [];
			end

			%% --------------------------------------------------------------------
			% SOLVE
			% ---------------------------------------------------------------------
			HJB = obj.solve_HJB(hours_bc_HJB);
			[KFE, Au] = obj.solve_KFE(HJB, hours_bc_KFE);
		    
			%% --------------------------------------------------------------------
			% COMPUTE WEALTH
			% ---------------------------------------------------------------------
			wealth = compute_wealth(KFE.g, obj.grids_KFE);
			fprintf('    --- MEAN WEALTH = %f ---\n\n', wealth)
		end
	end

	methods (Access=protected)
		function HJB = solve_HJB(obj, hours_bc)
			import HACTLib.computation.make_initial_guess

			[Vn, gguess] = make_initial_guess(obj.p, obj.grids_HJB,...
					obj.grids_KFE, obj.income);

			fprintf('    --- Iterating over HJB ---\n')
		    dst = 1e5;
			for nn	= 1:obj.p.HJB_maxiters
				[HJB, V_deriv_risky_asset_nodrift] = solver.find_policies(...
					obj.p, obj.income, obj.grids_HJB, Vn, hours_bc);

			    % Construct transition matrix 
		        [A, stationary] = obj.A_constructor_HJB.construct(HJB, Vn);

		        if obj.p.SDU
			    	risk_adj = {compute_risk_adjustment_for_nodrift_case(...
			    		obj.p, obj.grids_HJB, V_deriv_risky_asset_nodrift,...
			    		stationary, Vn)};
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

		    if (nn >= obj.p.HJB_maxiters)
		        error("HJB didn't converge");
		    end
		    
		    % Store value function and policies on both grids
		    HJB.Vn = Vn;
		end

		%% --------------------------------------------------------------------
	    % SOLVE KFE
		% ---------------------------------------------------------------------
		function [KFE, Au] = solve_KFE(obj, HJB, hours_bc)
			interp_decision = HACTLib.aux.interp_2d(obj.grids_KFE.b.vec,...
				obj.grids_KFE.a.vec, obj.grids_HJB.b.vec, obj.grids_HJB.a.vec);
			interp_decision = kron(speye(obj.income.ny*obj.p.nz), interp_decision);

			Vn = reshape(interp_decision * HJB.Vn(:),...
		    	[obj.p.nb_KFE, obj.p.na_KFE, obj.p.nz, obj.income.ny]);
		    KFE = solver.find_policies(obj.p, obj.income, obj.grids_KFE, Vn, hours_bc);
		    KFE.Vn = Vn;

			import HACTLib.computation.TransitionMatrixConstructor

			% True if returns should be treated as risky in the KFE
			returns_risk = (obj.p.sigma_r > 0) && (obj.p.retrisk_KFE == 1);
		    A_constructor_kfe = TransitionMatrixConstructor(obj.p,...
		    	obj.income, obj.grids_KFE, returns_risk);
		    Au = A_constructor_kfe.construct(KFE, KFE.Vn);
		    
		    
		    KFE.g = obj.kfe_solver.solve(Au);
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

function wealth = compute_wealth(g, grdKFE)
	iwealth= (g(:) .* grdKFE.trapezoidal.matrix(:))' * grdKFE.a.matrix(:);
	lwealth = (g(:) .* grdKFE.trapezoidal.matrix(:))' * grdKFE.b.matrix(:);
	wealth = iwealth + lwealth;
end