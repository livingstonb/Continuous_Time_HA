function [AYdiff,HJB,KFE,Au,grd,grdKFE] = solver(runopts,p,income,grd,grdKFE)
	% Solve for the steady state
    %
	% AYdiff is mean assets/income - target
    % HJB, KFE hold policy functions on the HJB & KFE grids
    % Au is transition matrix on KFE grid
    % grd, grdKFE are grids on HJB & KFE grids

    % keep track of number of mean asset iterations
	persistent iterAY
	if strcmp(runopts.RunMode,'Iterate') && isempty(iterAY)
		iterAY = 1;
	elseif strcmp(runopts.RunMode,'Iterate') && ~isempty(iterAY)
		iterAY = iterAY + 1;
	end

	fprintf('Solving for rho = %7.7f\n',p.rho)

	% ------------- unpack objects-----------------
	% asset grids
	nb = p.nb;
	nb_KFE = p.nb_KFE;
    nc = p.nc;
    nc_KFE = p.nc_KFE;

	% income
	ydist = income.ydist;
	y_mat = income.y.matrix;
	ny = income.ny;

	nz = p.nz;
	% ----------------------------------------------


	%% -----------------------------------------------------------
	% INITIALIZATION FOR HJB
	% ------------------------------------------------------------
	if numel(p.rhos) > 1
		rho_mat = reshape(p.rhos,[1 1 numel(p.rhos)]);
	else
		rho_mat =  p.rho;
	end

	% In some cases, replace value function guess and distribution with 
    % previous V and g`
    if ~ismember(runopts.RunMode,{'Iterate','Final'}) || (runopts.IterateRho==0)
        % First iteration, or not iterating, so don't load V_0
	    ConstructGuess = true;
    elseif strcmp(runopts.RunMode,'Iterate') && (iterAY==1)
        % Load results from rho lower bound outcome
        S = load([runopts.temp,'/V_lb_run',runopts.suffix,'.mat']);
        V_0 = S.V_0;
        gg0 = S.g(:);
        ConstructGuess = false;
    elseif strcmp(runopts.RunMode,'Iterate') && (iterAY==2)
    	% Load results from rho upper bound outcome
        S = load([runopts.temp,'/V_ub_run',runopts.suffix,'.mat']);
        V_0 = S.V_0;
        gg0 = S.g(:);
        ConstructGuess = false;
    else
		% Running in 'Final' mode after iterating or iterations are greater than 2
		% Load results from last fzero iteration
		S = load([runopts.temp,'/V0run',runopts.suffix,'.mat']);
		V_0 = S.V_0;
		gg0 = S.g(:);
		ConstructGuess = false;
    end

    if ConstructGuess
		% Initial guess
		[V_0, c_0] = solver.con_effort.find_guess(p,income,grd);
	    
		% adjust guess
		VbF_0 = zeros(nb,nz,ny);
		VbF_0(1:nb-1,:,:) = (V_0(2:nb,:,:) - V_0(1:nb-1,:,:)) ./ squeeze(grd.b.dF(1:nb-1,1));
		VbF_0 = reshape(VbF_0,[nb 1 nz ny]);
		VbF_0 = repmat(VbF_0,[1 nc 1 1]);

		VbB_0 = zeros(nb,nz,ny);
		VbB_0(2:nb,:,:) = (V_0(2:nb,:,:) - V_0(1:nb-1,:,:)) ./ squeeze(grd.b.dB(2:nb,1));
		VbB_0 = reshape(VbB_0,[nb 1 nz ny]);
		VbB_0 = repmat(VbB_0,[1 nc 1]);
	    
	    VbC_0 = zeros(nb,nz,ny);
	    bdC = zeros(nb,nz,ny);
	    bdC(2:nb-1,:,:) = repmat(grd.b.vec(3:nb) - grd.b.vec(1:nb-2),[1 ny]);
	    bdC(1,:,:) = bdC(2,:,:);
	    bdC(end,:,:) = bdC(end-1,:,:);
	    VbC_0(2:nb-1,:,:) = (V_0(3:nb,:,:) - V_0(1:nb-2,:,:)) ./ bdC(2:nb-1,:,:);
	    VbC_0 = reshape(VbC_0,[nb 1 nz ny]);
	    VbC_0 = repmat(VbC_0,[1 nc 1 1]);
	    VbC_0(nb,:,:,:) = VbB_0(nb,:,:,:);
	    VbC_0(1,:,:,:) = VbF_0(1,:,:,:);

	    bdot_0 = (p.r_b+p.deathrate * p.perfectannuities) * grd.b.matrix ...
	    	+ (1-p.wagetax) * y_mat - repmat(reshape(c_0,[nb 1 nz ny]),[1 nc 1 1]);
	    Vb_0 = VbF_0 .* (bdot_0>0)  + VbB_0 .* (bdot_0<0) + VbC_0 .* (bdot_0==0);

		bdot =  (p.r_b+p.deathrate*p.perfectannuities) * grd.b.matrix + (1-p.wagetax)*y_mat - grd.c.matrix;
	    bdot(1,:,:,:) = max(bdot(1,:,:,:),0);
	    
	    c = grd.c.matrix;
	    u = aux.u_fn(c,p.riskaver) - aux.con_effort.penalty(grd.b.matrix,p.penalty1,p.penalty2);
		V_0 = (1/(rho_mat+p.deathrate)) * (u + Vb_0 .* bdot);
	    
		% Initial distribution
	    gg0 = ones(nb_KFE,nc_KFE,nz,ny);
	    gg0 = gg0 .* permute(repmat(ydist,[1 nb_KFE nc_KFE nz]),[2 3 4 1]);
		gg0 = gg0 / sum(gg0(:));
		gg0 = gg0 ./ grdKFE.trapezoidal.matrix;
		gg0 = gg0(:);
	end
    
    % Create interpolation matrix between KFE grid and regular grid
	interp_decision = aux.interpTwoD(grdKFE.b.vec,grdKFE.c.vec,grd.b.vec,grd.c.vec);
	interp_decision = kron(speye(ny*nz),interp_decision);

	%% --------------------------------------------------------------------
    % FIND VALUE FUNCTION
	% ---------------------------------------------------------------------
	% Check if maximum number of iterations has been exceeded
	if iterAY > p.maxit_AY
		msgID = 'RHOITERATION:MaxIterExceeded';
	    msg = 'RHOITERATION:MaxIterExceeded';
	    iterException = MException(msgID,msg);
	    throw(iterException)
    end

    Vn = V_0;

    if numel(p.rhos) > 1
        rhocol = kron(p.rhos',ones(nb*na,1));
        rho_diag = spdiags(rhocol,0,nb*na*nz,nb*na*nz);
    else
        rho_diag = p.rho * speye(nb*na*nz);
	end
    
    fprintf('    --- Iterating over HJB ---\n')
    dst = 1e5;
    HJB.c = grd.c.matrix;
    HJB.s = (p.r_b+p.deathrate*p.perfectannuities) * grd.b.matrix + (1-p.wagetax) * y_mat - HJB.c;
	for nn	= 1:p.maxit_HJB
        % drift in c
        HJB.h = solver.con_effort.find_policies(p,income,grd,Vn);
        util = aux.u_fn(HJB.c,p.riskaver) - aux.con_effort.con_adj_cost(p,HJB.h)...
                - aux.con_effort.penalty(grd.b.matrix,p.penalty1,p.penalty2);

        A = solver.con_effort.construct_trans_matrix(p,income,grd,HJB); % A matrix

        Vn1 = solver.con_effort.update_value_function(p,income,A,util,HJB,Vn,rho_diag);
            
        dst = max(abs(Vn1(:)-Vn(:)));
        Vn = Vn1;

        if dst < p.crit_HJB
        	fprintf('\tHJB converged after %i iterations\n',nn);
		 	break
		elseif dst>10 && nn>1000
		 	% Not going to converge
		 	msgID = 'HJB:NotConverging';
		    msg = 'HJB:NotConverging';
		    HJBException = MException(msgID,msg);
		    throw(HJBException)
        elseif (nn==1) || (mod(nn,25) == 0)
            fprintf('\tHJB iteration = %i, distance = %e\n',nn,dst);
        end
        
    end
    
    if nn == p.maxit_HJB
        msgID = 'HJB:NoConvergence';
        msg = 'HJB:NoConvergence';
        HJBException = MException(msgID,msg);
        throw(HJBException)
    end
    
    HJB.Vn = Vn;

	%% --------------------------------------------------------------------
    % SOLVE KFE
	% ---------------------------------------------------------------------
	% put value functions onto KFE grid

	if (p.ny>1) && (p.nz>1)
		interp_grids = {grd.b.vec,grd.c.vec,grd.z.vec,income.y.vec};
	elseif p.ny > 1
		interp_grids = {grd.b.vec,grd.c.vec,income.y.vec};
	elseif p.nz > 1
		interp_grids = {grd.b.vec,grd.c.vec,grd.z.vec};
	else
		interp_grids = {grd.b.vec,grd.c.vec};
	end
	vinterp = griddedInterpolant(interp_grids,squeeze(HJB.Vn),'linear');

	if (p.ny>1) && (p.nz>1)
		v = vinterp(grdKFE.b.matrix(:),grdKFE.c.matrix(:),grd.z.matrix(:),income.y.matrixKFE(:));
	elseif p.ny > 1
		v = vinterp(grdKFE.b.matrix(:),grdKFE.c.matrix(:),income.y.matrixKFE(:));
	elseif p.nz > 1
		v = vinterp(grdKFE.b.matrix(:),grdKFE.c.matrix(:),grd.z.matrix(:));
	else
		v = vinterp(grdKFE.b.matrix(:),grdKFE.c.matrix(:));
	end

	KFE.Vn = reshape(v,[nb_KFE,nc_KFE,nz,income.ny]);

	KFE.h = solver.con_effort.find_policies(p,income,grdKFE,KFE.Vn);

    KFE.c = grdKFE.c.matrix;
     KFE.b = grdKFE.b.matrix;
    KFE.s = (p.r_b+p.deathrate*p.perfectannuities) * grdKFE.b.matrix...
            + (1-p.wagetax)*income.y.matrixKFE - grdKFE.c.matrix;

    % Construct A matrix for KFE
    Au = solver.con_effort.construct_trans_matrix(p,income,grdKFE,KFE);
    
    % solve
	KFE.g = solver.solveKFE(p,income,grdKFE,gg0,Au,'c');
    KFE.u = aux.u_fn(KFE.c,p.riskaver) - aux.con_effort.con_adj_cost(p,KFE.h)...
                - aux.con_effort.penalty(grdKFE.b.matrix,p.penalty1,p.penalty2);

    %% --------------------------------------------------------------------
	% STORE VALUE FUNCTION AND DISTRIBUTION FOR FUTURE ITERATIONS
	% ---------------------------------------------------------------------
	V_0 = Vn;
	g = KFE.g;
    if strcmp(runopts.RunMode,'Find rho_ub')
        save([runopts.temp,'V_ub_run',runopts.suffix,'.mat'],'V_0','g');
    elseif strcmp(runopts.RunMode,'Find rho_lb')
        save([runopts.temp,'V_lb_run',runopts.suffix,'.mat'],'V_0','g');
    elseif runopts.IterateRho == 1
    	% Save results for next fzero iteration
		save([runopts.temp,'V0run',runopts.suffix,'.mat'],'V_0','g')
	end

	%% --------------------------------------------------------------------
	% COMPUTE MEAN WEALTH
	% ---------------------------------------------------------------------
	Ewealth = (KFE.g(:) .* grdKFE.trapezoidal.matrix(:))' * grdKFE.b.matrix(:);

	AYdiff = Ewealth - p.targetAY;
	fprintf('    --- MEAN WEALTH = %f ---\n\n',Ewealth)
end