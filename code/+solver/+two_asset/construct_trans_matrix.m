function A = construct_trans_matrix(p, income, grids, model, modeltype, Vdiff_SDU, V)
    % Constructs transition matrix for the HJB when modeltype='HJB' and for
    % the KFE when modeltype='KFE'.
	na = numel(grids.a.vec);
	nb = numel(grids.b.vec);
    
    % policy functions
    s = model.s;
    d = model.d;

    if strcmp(modeltype,'KFE')
        y_mat = income.y.matrixKFE;
    elseif strcmp(modeltype,'HJB')
        y_mat = income.y.matrix;
    end
	ny = numel(income.y.vec);
    nz = p.nz;

    dim = nb * na * nz * ny;

	%% ----------------------------------------------
    % INITIALIZE
	% -----------------------------------------------
    chi = NaN(nb,na,nz,ny);
    yy = NaN(nb,na,nz,ny);
    zetta = NaN(nb,na,nz,ny);

    X = NaN(nb,na,nz,ny);
    Y = NaN(nb,na,nz,ny);
    Z = NaN(nb,na,nz,ny);
	
	%% --------------------------------------------------------------------
    % COMPUTE ASSET DRIFTS
	% ---------------------------------------------------------------------
    if strcmp(modeltype,'KFE')
        adriftB = min(d + grids.a.matrix * (p.r_a + p.deathrate*p.perfectannuities)...
                                                 + p.directdeposit * y_mat,0);
        adriftF = max(d + grids.a.matrix * (p.r_a + p.deathrate*p.perfectannuities)...
                                                 + p.directdeposit * y_mat,0);

        bdriftB = min(s - d - aux.two_asset.adj_cost_fn(d,grids.a.matrix,p),0);
        bdriftF = max(s - d - aux.two_asset.adj_cost_fn(d,grids.a.matrix,p),0);
    elseif strcmp(modeltype,'HJB')
        adriftB = min(d,0) + min(grids.a.matrix * (p.r_a + p.deathrate*p.perfectannuities) + p.directdeposit * y_mat,0);
        adriftF = max(d,0) + max(grids.a.matrix * (p.r_a + p.deathrate*p.perfectannuities) + p.directdeposit * y_mat,0);

        bdriftB = min(-d - aux.two_asset.adj_cost_fn(d,grids.a.matrix,p),0) + min(s,0);
        bdriftF = max(-d - aux.two_asset.adj_cost_fn(d,grids.a.matrix,p),0) + max(s,0);    
    end

	%% --------------------------------------------------------------------
    % ILLIQUID ASSET TRANSITIONS
	% ---------------------------------------------------------------------
    chi(:,2:na,:,:) = -adriftB(:,2:na,:,:)./grids.a.dB(:,2:na); 
    chi(:,1,:,:) = zeros(nb,1,nz,ny);
    yy(:,2:na-1,:,:) = adriftB(:,2:na-1,:,:)./grids.a.dB(:,2:na-1) ...
        - adriftF(:,2:na-1,:,:)./grids.a.dF(:,2:na-1); 
    yy(:,1,:,:) = - adriftF(:,1,:,:)./grids.a.dF(:,1); 
    yy(:,na,:,:) = adriftB(:,na,:,:)./grids.a.dB(:,na);
    zetta(:,1:na-1,:,:) = adriftF(:,1:na-1,:,:)./grids.a.dF(:,1:na-1); 
    zetta(:,na,:,:) = zeros(nb,1,nz,ny);

    centdiag = reshape(yy,nb*na,nz,ny);
    lowdiag = reshape(chi,nb*na,nz,ny);
    lowdiag = circshift(lowdiag,-nb);
    updiag = reshape(zetta,nb*na,nz,ny);
    updiag = circshift(updiag,nb);  

    centdiag = reshape(centdiag,dim,1);
    updiag   = reshape(updiag,dim,1);
    lowdiag  = reshape(lowdiag,dim,1);

    A = spdiags(centdiag,0,dim,dim) ...
        + spdiags(updiag,nb,dim,dim)...
        + spdiags(lowdiag,-nb,dim,dim);

	%% --------------------------------------------------------------------
    % LIQUID ASSET TRANSITIONS
	% ---------------------------------------------------------------------
    X(2:nb,:,:,:) = - bdriftB(2:nb,:,:,:)./grids.b.dB(2:nb,:); 
    X(1,:,:,:) = zeros(1,na,nz,ny);
    Y(2:nb-1,:,:,:) = bdriftB(2:nb-1,:,:,:)./grids.b.dB(2:nb-1,:) ...
        - bdriftF(2:nb-1,:,:,:)./grids.b.dF(2:nb-1,:); 
    Y(1,:,:,:) = - bdriftF(1,:,:,:)./grids.b.dF(1,:); 
    Y(nb,:,:,:) = bdriftB(nb,:,:,:)./grids.b.dB(nb,:);
    Z(1:nb-1,:,:,:) = bdriftF(1:nb-1,:,:,:)./grids.b.dF(1:nb-1,:); 
    Z(nb,:,:,:) = zeros(1,na,nz,ny);

    centdiag = reshape(Y,nb*na,nz,ny);
    lowdiag = reshape(X,nb*na,nz,ny);
    lowdiag = circshift(lowdiag,-1);
    updiag = reshape(Z,nb*na,nz,ny);
    updiag = circshift(updiag,1);

    centdiag = reshape(centdiag,dim,1);
    updiag   = reshape(updiag,dim,1);
    lowdiag  = reshape(lowdiag,dim,1);

    A = A + spdiags(centdiag,0,dim,dim) ...
        + spdiags(updiag,1,dim,dim)...
        + spdiags(lowdiag,-1,dim,dim);

    %% --------------------------------------------------------------------
    % RATE OF RETURN RISK: ONE ASSET
    % ---------------------------------------------------------------------
    if (p.sigma_r > 0) && (p.OneAsset == 1) && (strcmp(modeltype,'HJB') || (strcmp(modeltype,'KFE') && (p.retrisk_KFE==1)))
        term1 = ((grids.b.matrix * p.sigma_r) .^2) ./ (grids.b.dB + grids.b.dF);
        lowdiag = term1 ./ grids.b.dB;
        centdiag = - term1 .* (1 ./ grids.b.dF + 1 ./ grids.b.dB);
        updiag = term1 ./ grids.b.dF;

        % top of the grid
        term1 = (1/2) * (grids.b.matrix(nb,:,:,:) * p.sigma_r) .^ 2;
        lowdiag(nb,:,:,:) = term1 ./ (grids.b.dB(nb,:)).^2;
        centdiag(nb,:,:,:) = - term1 ./ (grids.b.dB(nb,:)).^2;
        updiag(nb,:,:,:) = 0;

        lowdiag = lowdiag(:);
        lowdiag = [lowdiag(2:end); 0];
        centdiag = centdiag(:);
        updiag = [0; updiag(1:end-1)];

        if p.SDU == 1
            % compute second difference, Vbb
            Vdiff2 = zeros(nb, na, nz, ny);
            grids_delta = grids.b.dF(:,1) + grids.b.dB(:,1);
            V_i_plus_1_term = 2 * V(3:nb,:,:,:) ./ (grids.b.dF(3:nb,:) .* grids_delta(3:nb));
            V_i_term = - 2 * V(2:nb-1,:,:,:) .* (1./grids.b.dF(2:nb-1,:) + 1./grids.b.dB(2:nb-1,:))...
                ./ grids_delta(2:nb-1);
            V_i_minus_1_term = 2 * V(1:nb-2,:,:,:) ./ (grids.b.dB(3:nb,:) .* grids_delta(1:nb-2));
            Vdiff2(2:nb-1,:,:,:) = V_i_plus_1_term + V_i_term + V_i_minus_1_term;
            Vdiff2(nb,:,:,:) = (V(nb-1,:,:,:) - V(nb,:,:,:)) ./ (grids.b.dB(nb,:) .^ 2);

            % compute adjustment term
            sdu_adj = 1 + ((p.invies-p.riskaver) / (1-p.invies)) * (Vdiff_SDU .^2)...
                ./ (V .* Vdiff2);
            sdu_adj = sdu_adj(:);
        
            lowdiag = sdu_adj .* lowdiag;
            centdiag = sdu_adj .* centdiag;
            updiag = sdu_adj .* updiag;
        end

        A = A + spdiags(centdiag,0,dim,dim);
        A = A + spdiags(updiag,1,dim,dim);
        A = A + spdiags(lowdiag,-1,dim,dim);
    elseif (p.sigma_r > 0) && (p.OneAsset == 0) && (strcmp(modeltype,'HJB') || (strcmp(modeltype,'KFE') && (p.retrisk_KFE==1)))
        term1 = ((grids.a.matrix * p.sigma_r) .^2) ./ (grids.a.dB + grids.a.dF);
        lowdiag = term1 ./ grids.a.dB;
        centdiag = - term1 .* (1 ./ grids.a.dF + 1 ./ grids.a.dB);
        updiag = term1 ./ grids.a.dF;
        
        % top of the grid
        term1 = (1/2) * (grids.a.matrix(:,na,:,:) * p.sigma_r) .^ 2;
        lowdiag(:,na,:,:) = term1 ./ (grids.a.dB(:,na)).^2;
        centdiag(:,na,:,:) = - term1 ./ (grids.a.dB(:,na)).^2;
        updiag(:,na,:,:) = 0;

        lowdiag = lowdiag(:);
        lowdiag = [lowdiag(nb+1:end); zeros(nb,1)];
        centdiag = centdiag(:);
        updiag = [zeros(nb,1); updiag(1:end-nb)];

        if p.SDU == 1
            % compute second difference, Vbb
            Vdiff2 = zeros(nb, na, nz, ny);
            grids_delta = grids.a.dF(:,1) + grids.a.dB(:,1);
            V_i_plus_1_term = 2 * V(3:na,:,:,:) ./ (grids.a.dF(3:na,:) .* grids_delta(3:na));
            V_i_term = - 2 * V(2:na-1,:,:,:) .* (1./grids.a.dF(2:na-1,:) + 1./grids.a.dB(2:na-1,:))...
                ./ grids_delta(2:na-1);
            V_i_minus_1_term = 2 * V(1:na-2,:,:,:) ./ (grids.a.dB(3:na,:) .* grids_delta(1:na-2));
            Vdiff2(2:na-1,:,:,:) = V_i_plus_1_term + V_i_term + V_i_minus_1_term;
            Vdiff2(na,:,:,:) = (V(na-1,:,:,:) - V(na,:,:,:)) ./ (grids.a.dB(na,:) .^ 2);

            % compute adjustment term
            sdu_adj = 1 + ((p.invies-p.riskaver) / (1-p.invies)) * (Vdiff_SDU .^2)...
                ./ (V .* Vdiff2);

            sdu_adj = sdu_adj(:);
        
            lowdiag = sdu_adj .* lowdiag;
            centdiag = sdu_adj .* centdiag;
            updiag = sdu_adj .* updiag;
        end

        A = A + spdiags(centdiag,0,dim,dim);
        A = A + spdiags(updiag,nb,dim,dim);
        A = A + spdiags(lowdiag,-nb,dim,dim);
    end

end