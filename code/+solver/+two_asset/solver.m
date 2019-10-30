function [AYdiff,HJB,KFE,Au] = solver(runopts,p,income,grd,grdKFE)
	% Solve for the steady state via the HJB and KFE equations, then
	% compute some relevant statistics.

    % keep track of number of mean asset iterations
	persistent iterAY
	if strcmp(runopts.RunMode,'Iterate') && isempty(iterAY)
		iterAY = 1;
	elseif strcmp(runopts.RunMode,'Iterate') && ~isempty(iterAY)
		iterAY = iterAY + 1;
	end

	fprintf('Solving for rho = %7.7f\n',p.rho)

    na = p.na;
    nb = p.nb;
	na_KFE = p.na_KFE;
	nb_KFE = p.nb_KFE;
	ny = numel(income.y.vec);
	nz = p.nz;

	% Create interpolation matrix between KFE grid and regular grid
	interp_decision = aux.interpTwoD(grdKFE.b.vec,grdKFE.a.vec,grd.b.vec,grd.a.vec);
	interp_decision = kron(speye(ny*nz),interp_decision);

	% Returns grid
	r_b_mat = p.r_b .* (grd.b.matrix>=0) +  p.r_b_borr .* (grd.b.matrix<0);

	%% --------------------------------------------------------------------
	% INITIALIZATION FOR HJB
	% ---------------------------------------------------------------------
	if numel(p.rhos) > 1
		rho_mat = reshape(p.rhos,[1 1 numel(p.rhos)]);
	else
		rho_mat =  p.rho;
	end

	% Initial guess
	c_0 = (1-p.directdeposit - p.wagetax) * income.y.matrix ...
            + (p.r_a + p.deathrate*p.perfectannuities) * grd.a.matrix...
            + (r_b_mat + p.deathrate*p.perfectannuities) .* grd.b.matrix + p.transfer;

    if p.SDU == 0
    	V_0	= 1 ./ (rho_mat + p.deathrate) .* aux.u_fn(c_0,p.riskaver);
    else
    	V_0 = aux.u_fn(c_0, p.invies);
    end
    
%     % Attempt at more accurate guess (seems to fail)
%     inc_trans = kron(income.ytrans,speye(p.nb*p.na));
%     inc_trans = repmat(inc_trans,nz);
%     u_0 = aux.u_fn(c_0(:),p.riskaver);
%     if numel(p.rhos) == 1
%         V_0 = (p.rho*speye(p.nb*p.na*ny*p.nz) - inc_trans) \ u_0;
%     else
%         V_0 = (kron(diag(p.rhos'),speye(p.nb*p.na*ny)) - inc_trans) \ u_0;
%     end

    V_0 = reshape(V_0,[nb na nz ny]);

	% Initial distribution
    gg0 = ones(nb_KFE,na_KFE,nz,ny);
    gg0 = gg0 .* permute(repmat(income.ydist,[1 nb_KFE na_KFE nz]),[2 3 4 1]);
    if p.OneAsset == 1
        gg0(:,grdKFE.a.vec>0,:,:) = 0;
    end
    gg0 = gg0 / sum(gg0(:));
	gg0 = gg0 ./ grdKFE.trapezoidal.matrix;
	gg = gg0;

	%% --------------------------------------------------------------------
    % SOLVE HJB
	% ---------------------------------------------------------------------
	% Check if maximum number of iterations has been exceeded
	if iterAY > p.maxit_AY
		msgID = 'RHOITERATION:MaxIterExceeded';
	    msg = 'RHOITERATION:MaxIterExceeded';
	    iterException = MException(msgID,msg);
	    throw(iterException)
    end

    % In some cases, replace value function guess and distribution with 
    % previous V and g`
    if strcmp(runopts.RunMode,'NoRisk')
        Vn = V_0;
    elseif ~ismember(runopts.RunMode,{'Iterate','Final'}) || (runopts.IterateRho==0)
        % First iteration, or not iterating, so don't load V_0
	    Vn	= V_0;	
    elseif strcmp(runopts.RunMode,'Iterate') && (iterAY==1)
        % Load results from rho lower bound outcome
        S = load([runopts.temp,'/V_lb_run',runopts.suffix,'.mat']);
        Vn = S.V_0;
        gg = S.g(:);
    elseif strcmp(runopts.RunMode,'Iterate') && (iterAY==2)
    	% Load results from rho upper bound outcome
        S = load([runopts.temp,'/V_ub_run',runopts.suffix,'.mat']);
        Vn = S.V_0;
        gg = S.g(:);
    else
		% Running in 'Final' mode after iterating or iterations are greater than 2
		% Load results from last fzero iteration
		S = load([runopts.temp,'/V0run',runopts.suffix,'.mat']);
		Vn = S.V_0;
		gg = S.g(:);
    end

    fprintf('    --- Iterating over HJB ---\n')
    dst = 1e5;
	for nn	= 1:p.maxit_HJB
	  	HJB = solver.two_asset.find_policies(p,income,grd,Vn);

	    % CONSTRUCT TRANSITION MATRIX 
        A = solver.two_asset.construct_trans_matrix(p,income,grd,HJB,'HJB');

	    % UPDATE VALUE FUNCTION
	    Vn1 = solver.two_asset.solveHJB(p,A,income,Vn,HJB.u,nn);

	    % Check for convergence
	    Vdiff = Vn1 - Vn;
	    Vn = Vn1;
	    dst = max(abs(Vdiff(:)));
        if (nn==1) || (mod(nn,25)==0)
	    	fprintf('\tHJB iteration = %i, distance = %e\n',nn,dst);
        end

        if dst<p.crit_HJB
	        fprintf('\tHJB converged after %i iterations\n',nn);
		 	break
		elseif dst>10 && nn>500
		 	% Not going to converge, throw exception
		 	msgID = 'HJB:NotConverging';
		    msg = 'HJB:NotConverging';
		    HJBException = MException(msgID,msg);
		    throw(HJBException)
		end
    end
    
    if (nn >= p.maxit_HJB)
        error("HJB didn't converge");
    end
    
    % Store value function and polices on both grids
    HJB.Vn = Vn;
    KFE_Vn = reshape(interp_decision*Vn(:),nb_KFE,na_KFE,nz,ny);
    KFE = solver.two_asset.find_policies(p,income,grdKFE,KFE_Vn);
    KFE.Vn = KFE_Vn;

	%% --------------------------------------------------------------------
    % SOLVE KFE
	% ---------------------------------------------------------------------
    Au = solver.two_asset.construct_trans_matrix(p,income,grdKFE,KFE,'KFE');

    dim2Identity = 'a';
	g = solver.solveKFE(p,income,grdKFE,gg,Au,dim2Identity);

	%% --------------------------------------------------------------------
	% STORE VALUE FUNCTION AND DISTRIBUTION FOR FUTURE ITERATIONS
	% ---------------------------------------------------------------------
	V_0 = Vn;
    if strcmp(runopts.RunMode,'Find rho_ub')
        save([runopts.temp,'V_ub_run',runopts.suffix,'.mat'],'V_0','g');
    elseif strcmp(runopts.RunMode,'Find rho_lb')
        save([runopts.temp,'V_lb_run',runopts.suffix,'.mat'],'V_0','g');
    elseif runopts.IterateRho == 1
    	% Save results for next fzero iteration
		save([runopts.temp,'V0run',runopts.suffix,'.mat'],'V_0','g')
	end

	%% --------------------------------------------------------------------
	% COMPUTE WEALTH
	% ---------------------------------------------------------------------
	Eillwealth = (g(:) .* grdKFE.trapezoidal.matrix(:))' * grdKFE.a.matrix(:);
	Eliqwealth = (g(:) .* grdKFE.trapezoidal.matrix(:))' * grdKFE.b.matrix(:);
	Etotwealth = Eillwealth + Eliqwealth;
    
    KFE.g = g;

	AYdiff = Etotwealth - p.targetAY;
	fprintf('    --- MEAN WEALTH = %f ---\n\n',Etotwealth)
    
end
