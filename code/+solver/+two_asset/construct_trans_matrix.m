function A = construct_trans_matrix(p,income,grids,model,modeltype)
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

    centdiag = reshape(centdiag,nb*na*ny*nz,1);
    updiag   = reshape(updiag,nb*na*ny*nz,1);
    lowdiag  = reshape(lowdiag,nb*na*ny*nz,1);

    A = spdiags(centdiag,0,nb*na*ny*nz,nb*na*ny*nz) ...
        + spdiags(updiag,nb,nb*na*ny*nz,nb*na*ny*nz)...
        + spdiags(lowdiag,-nb,nb*na*ny*nz,nb*na*ny*nz);

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

    centdiag = reshape(centdiag,nb*na*ny*nz,1);
    updiag   = reshape(updiag,nb*na*ny*nz,1);
    lowdiag  = reshape(lowdiag,nb*na*ny*nz,1);

    A = A + spdiags(centdiag,0,nb*na*ny*nz,nb*na*ny*nz) ...
        + spdiags(updiag,1,nb*na*ny*nz,nb*na*ny*nz)...
        + spdiags(lowdiag,-1,nb*na*ny*nz,nb*na*ny*nz);
end