function [HJB, KFE, Au] = solver(runopts, p, income, grd, grdKFE)
	% Solve for the steady state via the HJB and KFE equations, then
	% compute some relevant statistics.

	% Parameters
	% ----------
	% runopts : a structure of boolean variables indicating various
	%	options
	%
	% p : a Params object
	%
	% income : an Income object
	%
	% grd : a Grid object for the HJB
	%
	% grdKFE : a Grid object for the KFE
	%
	% Returns
	% -------
	% HJB : a structure containing policy functions and the value
	%	function on the HJB grid, shape (nb, na, nz, ny)
	%
	% KFE : a structure containing policy functions and the value
	%	function on the KFE grid, shape (nb_KFE, na_KFE, nz, ny)
	%
	% Au : the transition matrix for the KFE, shape
	%	(nb_KFE*na_KFE*nz*ny, nb_KFE*na_KFE*nz*ny)

    % keep track of number of mean asset iterations
	persistent iterAY
	if strcmp(runopts.RunMode,'Iterate') && isempty(iterAY)
		iterAY = 1;
	elseif strcmp(runopts.RunMode,'Iterate') && ~isempty(iterAY)
		iterAY = iterAY + 1;
	end

	fprintf('Solving for rho = %7.7f\n',p.rho)

	%% --------------------------------------------------------------------
	% INITIALIZATION
	% ---------------------------------------------------------------------
    na = p.na;
    nb = p.nb;
	na_KFE = p.na_KFE;
	nb_KFE = p.nb_KFE;
	ny = numel(income.y.vec);
	nz = p.nz;

	% create interpolation matrix between KFE grid and regular grid
	interp_decision = aux.interpTwoD(grdKFE.b.vec,grdKFE.a.vec,grd.b.vec,grd.a.vec);
	interp_decision = kron(speye(ny*nz),interp_decision);
	
	% initial guess for value function and distribution
    [Vn, gg] = solver.make_initial_guess(p, grd, grdKFE, income);
    
	%% --------------------------------------------------------------------
    % SOLVE HJB
	% ---------------------------------------------------------------------
	check_if_max_iters_exceeded(iterAY, p);

    returns_risk = (p.sigma_r > 0);
    A_Constructor = solver.A_Matrix_Constructor(p, income, grd, 'HJB', returns_risk);

    fprintf('    --- Iterating over HJB ---\n')
    dst = 1e5;
	for nn	= 1:p.maxit_HJB
	  	[HJB, V_deriv_risky_asset_nodrift] = solver.find_policies(p,income,grd,Vn);

	    % construct transition matrix 
        [A, stationary] = A_Constructor.construct(HJB, Vn);

    	risk_adj = compute_risk_adjustment_for_nodrift_case(...
    		p, grd, V_deriv_risky_asset_nodrift, stationary, Vn);
        
		% update value function
		Vn1 = solver.solveHJB(p, A, income, Vn, HJB.u, nn, risk_adj);

	    % check for convergence
	    Vdiff = Vn1 - Vn;
	    Vn = Vn1;
	    dst = max(abs(Vdiff(:)));
        if (nn==1) || (mod(nn,25)==0)
	    	fprintf('\tHJB iteration = %i, distance = %e\n',nn,dst);
        end

        if dst < p.crit_HJB
	        fprintf('\tHJB converged after %i iterations\n',nn);
		 	break
		end
		check_if_not_converging(dst, nn);
    end
    
    if (nn >= p.maxit_HJB)
        error("HJB didn't converge");
    end
    
    % store value function and policies on both grids
    HJB.Vn = Vn;
    KFE_Vn = reshape(interp_decision*Vn(:),nb_KFE,na_KFE,nz,ny);
    KFE = solver.find_policies(p,income,grdKFE,KFE_Vn);
    KFE.Vn = KFE_Vn;

	%% --------------------------------------------------------------------
    % SOLVE KFE
	% ---------------------------------------------------------------------
	% true if returns should be treated as risky in the KFE
	returns_risk = (p.sigma_r > 0) && (p.retrisk_KFE == 1);
    A_Constructor_KFE = solver.A_Matrix_Constructor(p, income, grdKFE, 'KFE', returns_risk);
    Au = A_Constructor_KFE.construct(KFE, KFE.Vn);
    
    kfe_solver = solver.KFESolver(p, income, grdKFE);
    KFE.g = kfe_solver.solve(Au);
    
	%% --------------------------------------------------------------------
	% COMPUTE WEALTH
	% ---------------------------------------------------------------------
	wealth = compute_wealth(KFE.g, grdKFE);
	fprintf('    --- MEAN WEALTH = %f ---\n\n',wealth)
end

function check_if_max_iters_exceeded(iterAY, p)
	if iterAY > p.maxit_AY
		msgID = 'RHOITERATION:MaxIterExceeded';
	    msg = 'RHOITERATION:MaxIterExceeded';
	    iterException = MException(msgID,msg);
	    throw(iterException)
    end
end

function check_if_not_converging(dst, nn)
	if dst>10 && nn>500
	 	% Not going to converge, throw exception
	 	msgID = 'HJB:NotConverging';
	    msg = 'HJB:NotConverging';
	    HJBException = MException(msgID,msg);
	    throw(HJBException)
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

function wealth = compute_wealth(g, grdKFE)
	iwealth= (g(:) .* grdKFE.trapezoidal.matrix(:))' * grdKFE.a.matrix(:);
	lwealth = (g(:) .* grdKFE.trapezoidal.matrix(:))' * grdKFE.b.matrix(:);
	wealth = iwealth + lwealth;
end