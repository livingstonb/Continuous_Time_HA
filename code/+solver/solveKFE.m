function g = solveKFE(p,income,grdKFE,gg,A,na_KFEIdentity)
	% solveKFE() iterates over the Kolmogorov Forward
	% Equation to find the equilibrium distribution over states

	nz = p.nz;
	nb_KFE = p.nb_KFE;
	na_KFE = p.na_KFE;
	ny = numel(income.y.vec);

	if p.iterateKFE == 0
        gg = ones(nb_KFE*na_KFE*nz, 1) .* income.ydist(:)' / (nb_KFE*na_KFE*nz);
        gg = gg(:);
        
		inctrans = kron(sparse(income.ytrans), speye(nb_KFE*na_KFE*nz));
% 		Ap_extended = sparse([(A+inctrans)'; grdKFE.trapezoidal.matrix(:)']);
        Ap_extended = sparse([(A+inctrans)'; ones(1, nb_KFE*na_KFE*nz*ny)]);
		RHS = sparse([zeros(nb_KFE*na_KFE*nz*ny, 1); 1]);

		gg = Ap_extended \ RHS;
%         gg = lsqr(Ap_extended, RHS, 1e-12, 1000000, [], [], gg(:));
        gg = gg ./ grdKFE.trapezoidal.matrix(:);
        gg = full(gg);
	else
		gg = gg(:);

	    gg1_denom = cell(1,ny);
	    for iy = 1:ny
	        Ak = A(1+(iy-1)*(nb_KFE*na_KFE*nz):iy*(nb_KFE*na_KFE*nz),1+(iy-1)*(nb_KFE*na_KFE*nz):iy*(nb_KFE*na_KFE*nz));
	        gg1_denom{iy} = (speye(nb_KFE*na_KFE*nz) - p.delta_KFE * Ak'...
		   		- p.delta_KFE * (income.ytrans(iy,iy) - p.deathrate) * speye(nb_KFE*na_KFE*nz));
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

		    if strcmp(na_KFEIdentity,'c')
		   		gg_tilde_cz_k = reshape(gg_tilde,[nb_KFE na_KFE*nz ny]);
		    end
		    gg1 = NaN(nb_KFE*na_KFE*nz,ny);

		    for iy = 1:ny    
	            gk_sum = sum(repmat(ytrans0p(iy,:),nb_KFE*na_KFE*nz,1) ...
	            	.* reshape(gg_tilde,nb_KFE*na_KFE*nz,ny),2);

	            if (p.Bequests == 1) && (p.ResetIncomeUponDeath == 1)
	                deathg = p.deathrate * income.ydist(iy) * sum(reshape(gg_tilde,[],ny),2);
	            elseif (p.Bequests == 1) && (p.ResetIncomeUponDeath == 0)
	                deathg = p.deathrate * gg_tilde(1+(iy-1)*(nb_KFE*na_KFE*nz):iy*(nb_KFE*na_KFE*nz));
	            elseif (p.Bequests == 0) && (p.ResetIncomeUponDeath == 1)
	                deathg = sparse(nb_KFE*na_KFE*nz,1);
	                deathg(1:nb_KFE*na_KFE:end) = p.deathrate * income.ydist(iy) * (1/nz);
	            elseif (p.Bequests == 0) && (p.ResetIncomeUponDeath == 0)
	                deathg = sparse(nb_KFE*na_KFE*nz,1);
	                if strcmp(na_KFEIdentity,'a')
	                	deathg(grdKFE.loc0b0a:nb_KFE*na_KFE:end) = p.deathrate * income.ydist(iy) * (1/nz);
	                elseif strcmp(na_KFEIdentity,'c')
	                	deathg(grdKFE.loc0b:nb_KFE:end) = p.deathrate * sum(gg_tilde_cz_k(:,:,iy),1);
	                end
	            end

	            gg1(:,iy) = gg1_denom{iy}*(gg_tilde(1+(iy-1)*(nb_KFE*na_KFE*nz):iy*(nb_KFE*na_KFE*nz))...
	                                     + p.delta_KFE*gk_sum + p.delta_KFE*deathg);
	        end

		    gg1 = reshape(gg1,nb_KFE*na_KFE*nz*ny,1);
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
	end

	g = reshape(gg,nb_KFE,na_KFE,nz,ny);
end