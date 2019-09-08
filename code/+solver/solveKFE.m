function g = solveKFE(p,income,grdKFE,gg,A,dim2Identity)
	% solveKFE() iterates over the Kolmogorov Forward
	% Equation to find the equilibrium distribution over states
	%
	% dim2Identity is either 'a' for the two-asset model,
	% or 'c' for the model with consumption adjustment costs

	if strcmp(dim2Identity,'a')
		dim2 = p.na_KFE;
		ResetIncomeUponDeath = p.ResetIncomeUponDeath;
	elseif strcmp(dim2Identity,'c')
		dim2 = p.nc_KFE;
		ResetIncomeUponDeath = 0;
	else
		error('no second dimension in grids')
	end

	nz = p.nz;
	nb_KFE = p.nb_KFE;
	ny = numel(income.y.vec);

	gg = gg(:);

    gg1_denom = cell(1,ny);
    for iy = 1:ny
        Ak = A(1+(iy-1)*(nb_KFE*dim2*nz):iy*(nb_KFE*dim2*nz),1+(iy-1)*(nb_KFE*dim2*nz):iy*(nb_KFE*dim2*nz));
        gg1_denom{iy} = (speye(nb_KFE*dim2*nz) - p.delta_KFE * Ak'...
	   		- p.delta_KFE * (income.ytrans(iy,iy) - p.deathrate) * speye(nb_KFE*dim2*nz));
        gg1_denom{iy} = inverse(gg1_denom{iy});
    end

    % transition matrix with diagonal killed
	ytrans0  = income.ytrans - diag(diag(income.ytrans)); 
	ytrans0p = ytrans0';

	iter = 0;
	dst = 1e5;
	fprintf('    --- Iterating over KFE ---\n')
	while iter <= p.maxit_KFE && dst >= p.crit_KFE
		iter = iter + 1;

	    gg_tilde = grdKFE.trapezoidal.diagm * gg;

	    if strcmp(dim2Identity,'c')
	   		gg_tilde_cz_k = reshape(gg_tilde,[nb_KFE dim2*nz ny]);
	    end
	    gg1 = NaN(nb_KFE*dim2*nz,ny);

	    for iy = 1:ny    
            gk_sum = sum(repmat(ytrans0p(iy,:),nb_KFE*dim2*nz,1) ...
            	.* reshape(gg_tilde,nb_KFE*dim2*nz,ny),2);

            if (p.Bequests == 1) && (ResetIncomeUponDeath == 1)
                deathg = p.deathrate * income.ydist(iy) * sum(reshape(gg_tilde,[],ny),2);
            elseif (p.Bequests == 1) && (ResetIncomeUponDeath == 0)
                deathg = p.deathrate * gg_tilde(1+(iy-1)*(nb_KFE*dim2*nz):iy*(nb_KFE*dim2*nz));
            elseif (p.Bequests == 0) && (ResetIncomeUponDeath == 1)
                deathg = sparse(nb_KFE*dim2*nz,1);
                deathg(1:nb_KFE*dim2:end) = p.deathrate * income.ydist(iy) * (1/nz);
            elseif (p.Bequests == 0) && (ResetIncomeUponDeath == 0)
                deathg = sparse(nb_KFE*dim2*nz,1);
                if strcmp(dim2Identity,'a')
                	deathg(1:nb_KFE*dim2:end) = p.deathrate * income.ydist(iy) * (1/nz);
                elseif strcmp(dim2Identity,'c')
                	deathg(1:nb_KFE:end) = p.deathrate * sum(gg_tilde_cz_k(:,:,iy),1);
                end
            end

            gg1(:,iy) = gg1_denom{iy}*(gg_tilde(1+(iy-1)*(nb_KFE*dim2*nz):iy*(nb_KFE*dim2*nz))...
                                     + p.delta_KFE*gk_sum + p.delta_KFE*deathg);
        end

	    gg1 = reshape(gg1,nb_KFE*dim2*nz*ny,1);
	    gg1 = gg1 ./ sum(gg1);
	    gg1 = grdKFE.trapezoidal.diagm \ gg1;

        dst = max(abs(gg1(:) - gg(:)));

        if (dst>10000) && (iter>2000)
        	msgID = 'KFE:NotConverging';
		    msg = 'KFE:NotConverging';
		    KFEException = MException(msgID,msg);
		    throw(KFEException)
        end
        
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

	g = reshape(gg,nb_KFE,dim2,nz,ny);
end