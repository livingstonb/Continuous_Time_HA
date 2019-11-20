function g = solveKFE(p, income, grdKFE, gg, A, V)
	% iterates over the Kolmogorov Forward
	% Equation to find the equilibrium distribution over states

	% Parameters
	% ----------
	% p : a Params object
	%
	% income : an Income object
	%
	% grdKFE : the asset grids for the KFE
	%
	% gg : the initial equilibrium distribution of states
	%
	% A : KFE transition matrix, of shape
	%	(nb_KFE*na_KFE*nz*ny, nb_KFE*na_KFE*nz*ny)
	%
	% V : value function, shape (nb_KFE, na_KFE, nz, ny)
	%
	% Returns
	% -------
	% g : the equilibrium distribution, of shape
	%	(nb_KFE, na_KFE, nz, ny)

	nz = p.nz;
	nb_KFE = p.nb_KFE;
	na_KFE = p.na_KFE;
	ny = numel(income.y.vec);

	if (p.SDU == 1) && (ny > 1)
        ez_adj = income.SDU_income_risk_adjustment(p, V);
        ez_adj_no_diag = ez_adj;
        for k = 1:ny
        	ez_adj_no_diag(:,k,k) = 0;
        end
    else
    	ez_adj = [];
    end

	if p.iterateKFE == 0
        gg = ones(nb_KFE*na_KFE*nz, 1) .* income.ydist(:)' / (nb_KFE*na_KFE*nz);
        gg = gg(:);
   
        A = A + income.sparse_income_transitions(p, ez_adj, 'KFE');
        Ap_extended = sparse([(A+inctrans)'; ones(1, nb_KFE*na_KFE*nz*ny)]);
		RHS = sparse([zeros(nb_KFE*na_KFE*nz*ny, 1); 1]);

		gg = Ap_extended \ RHS;
%         gg = lsqr(Ap_extended, RHS, 1e-12, 1000000, [], [], gg(:));
        gg = gg ./ grdKFE.trapezoidal.matrix(:);
        gg = full(gg);
	else
		gg = gg(:);

		% compute the LHS of the KFE
		KFE_LHS = KFE_matrix_divisor(p, income, A, ez_adj);

	    % transition matrix with diagonal killed
		ytrans0  = income.ytrans - diag(diag(income.ytrans)); 
		ytrans0p = ytrans0';

		iter = 0;
		dst = 1e5;
		fprintf('    --- Iterating over KFE ---\n')
		while iter <= p.maxit_KFE && dst >= p.crit_KFE
			iter = iter + 1;

		    gg_tilde = grdKFE.trapezoidal.diagm * gg;

		    gg1 = NaN(nb_KFE*na_KFE*nz,ny);

		    for iy = 1:ny    
		    	if (p.SDU == 0) || (ny == 1)
		            gk_sum = sum(repmat(ytrans0p(iy,:),nb_KFE*na_KFE*nz,1) ...
		            	.* reshape(gg_tilde,nb_KFE*na_KFE*nz,ny),2);
		        else
		        	gk_sum = sum(ez_adj_no_diag(:,:,k) .* reshape(gg_tilde,nb_KFE*na_KFE*nz,ny),2);
		        end
		        	
	            if (p.Bequests == 1) && (p.ResetIncomeUponDeath == 1)
	                deathg = p.deathrate * income.ydist(iy) * sum(reshape(gg_tilde,[],ny),2);
	            elseif (p.Bequests == 1) && (p.ResetIncomeUponDeath == 0)
	                deathg = p.deathrate * gg_tilde(1+(iy-1)*(nb_KFE*na_KFE*nz):iy*(nb_KFE*na_KFE*nz));
	            elseif (p.Bequests == 0) && (p.ResetIncomeUponDeath == 1)
	                deathg = sparse(nb_KFE*na_KFE*nz,1);
	                deathg(1:nb_KFE*na_KFE:end) = p.deathrate * income.ydist(iy) * (1/nz);
	            elseif (p.Bequests == 0) && (p.ResetIncomeUponDeath == 0)
	                deathg = sparse(nb_KFE*na_KFE*nz,1);
	                deathg(grdKFE.loc0b0a:nb_KFE*na_KFE:end) = p.deathrate * income.ydist(iy) * (1/nz);
  	            end

	            gg1(:,iy) = KFE_LHS{iy}*(gg_tilde(1+(iy-1)*(nb_KFE*na_KFE*nz):iy*(nb_KFE*na_KFE*nz))...
	                                     + p.delta_KFE*gk_sum + p.delta_KFE*deathg);
	        end

		    gg1 = reshape(gg1,nb_KFE*na_KFE*nz*ny,1);
		    gg1 = gg1 ./ sum(gg1);
		    gg1 = grdKFE.trapezoidal.diagm \ gg1;

	        dst = max(abs(gg1(:) - gg(:)));
	        check_if_not_converging(dst, iter);
	        
		    if (iter==1) || (mod(iter,100) == 0)
		        fprintf('\tKFE iteration  = %i, distance = %e\n',iter,dst);
		    end

		    gg = gg1;
		end

		if dst < p.crit_KFE
		    fprintf('\tKFE converged after %i iterations\n',iter);
		elseif dst >= p.crit_KFE
			error('KFE did not converge')
		end
	end

	g = reshape(gg,nb_KFE,na_KFE,nz,ny);
end

function LHS = KFE_matrix_divisor(p, income, A, ez_adj)
	% computes the left-hand-side of the KFE

	% Parameters
	% ----------
	% p : a Params object
	%
	% income : an Income object
	%
	% A : the KFE transition matrix
	%
	% ez_adj : the risk-adjusted income transitions
	%
	% Returns
	% -------
	% LHS : a cell array of operators B_k s.t. B_k * RHS_k
	%	returns the k-th income section of the equilibrium distribution,
	%	which is LHS_k \ RHS_k

	LHS = cell(1, income.ny);
	for k = 1:income.ny
		i1 = 1 + (k-1) * (p.nb_KFE*p.na_KFE*p.nz);
		i2 = k * (p.nb_KFE*p.na_KFE*p.nz);
		Ak = A(i1:i2, i1:i2);

		if (p.SDU == 0) || (income.ny == 1)
			LHS{k} = (speye(p.nb_KFE*p.na_KFE*p.nz) - p.delta_KFE * Ak'...
	   			- p.delta_KFE * (income.ytrans(k,k) - p.deathrate) * speye(p.nb_KFE*p.na_KFE*p.nz));
		else
			LHS{k} = (1+p.delta_KFE*p.deathrate)*speye(p.nb_KFE*p.na_KFE*p.nz) - p.delta_KFE * Ak'...
	   			- p.delta_KFE * spdiags(ez_adj(:,k,k), 0, p.nb_KFE*p.na_KFE*p.nz, p.nb_KFE*p.na_KFE*p.nz);
		end
		LHS{k} = inverse(LHS{k});
	end
end

function check_if_not_converging(dst, iter)
	if (dst>10000) && (iter>2000)
    	msgID = 'KFE:NotConverging';
	    msg = 'KFE:NotConverging';
	    KFEException = MException(msgID,msg);
	    throw(KFEException)
    end
end