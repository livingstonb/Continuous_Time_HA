function [A, stationary] = construct_trans_matrix(p, income, grids, model, modeltype, Vdiff_SDU, V)
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

    stationary = [];
	
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

%     adriftB = max(adriftB, -100);
%     adriftF = min(adriftF, 100);
%     bdriftB = max(bdriftB, -100);
%     bdriftF = min(bdriftF, 100);

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

        lowdiag = circshift(reshape(lowdiag, nb*na, nz, ny), -nb);
        updiag = circshift(reshape(updiag, nb*na, nz, ny), nb);

        lowdiag = lowdiag(:);
        updiag = updiag(:);
        centdiag = centdiag(:);

        if p.SDU == 1
            % compute second difference, Vbb
            Vdiff2 = zeros(nb, na, nz, ny);
            dF_adj = grids.b.dF(:,1);
            dF_adj(nb) = grids.b.dB(nb,1);
            grids_delta = dF_adj + grids.b.dB(:,1);
            V_i_plus_1_term = 2 * V(3:nb,:,:,:) ./ (dF_adj(3:nb) .* grids_delta(3:nb));
            V_i_term = - 2 * V(2:nb-1,:,:,:) .* (1./dF_adj(2:nb-1) + 1./grids.b.dB(2:nb-1,:))...
                ./ grids_delta(2:nb-1);
            V_i_minus_1_term = 2 * V(1:nb-2,:,:,:) ./ (grids.b.dB(1:nb-2,:) .* grids_delta(1:nb-2));
            Vdiff2(2:nb-1,:,:,:) = V_i_plus_1_term + V_i_term + V_i_minus_1_term;
            Vdiff2(nb,:,:,:) = (V(nb-1,:,:,:) - V(nb,:,:,:)) ./ (grids.b.dB(nb,:) .^ 2);

            % compute adjustment term
            if p.invies == 1
                sdu_adj = 1 + (1-p.riskaver) * (Vdiff_SDU .^2) ./ Vdiff2;
            else
                sdu_adj = 1 + ((p.invies-p.riskaver) / (1-p.invies)) * (Vdiff_SDU .^2)...
                    ./ (V .* Vdiff2);
            end
            sdu_adj(1,:,:,:) = 0;
            sdu_adj = sdu_adj(:);
        
            lowdiag = sdu_adj .* lowdiag;
            centdiag = sdu_adj .* centdiag;
            updiag = sdu_adj .* updiag;
        end

        A = A + spdiags(centdiag,0,dim,dim);
        A = A + spdiags(updiag,1,dim,dim);
        A = A + spdiags(lowdiag,-1,dim,dim);
    elseif (p.sigma_r > 0) && (p.OneAsset == 0) && (strcmp(modeltype,'HJB') || (strcmp(modeltype,'KFE') && (p.retrisk_KFE==1)))
       
        % % 2nd derivative term, Vaa
        % deltas = grids.a.dB + grids.a.dF;
        % deltas(:, na) = 2 * grids.a.dB(:, na);

        % risk_term = (grids.a.matrix * p.sigma_r) .^2; % (1/2) cancels out
        % lowdiag = risk_term ./ (grids.a.dB .* deltas);

        % centdiag = - risk_term .* (1 ./ grids.a.dF + 1 ./ grids.a.dB) ./ deltas;
        % centdiag(:, na, :, :) = -risk_term(:, na, :, :) ./ ( grids.a.dB(:, na) .* deltas(:, na) );

        % updiag = risk_term ./ (grids.a.dF .* deltas);
        % updiag(:, na, :, :) = 0;

        % % lowdiag = circshift(reshape(lowdiag, nb*na, nz, ny), -nb);
        % % updiag = circshift(reshape(updiag, nb*na, nz, ny), nb);

        % % A = A + spdiags(centdiag(:),0,dim,dim);
        % % A = A + spdiags(updiag(:),nb,dim,dim);
        % % A = A + spdiags(lowdiag(:),-nb,dim,dim);
        
        % %% THIS SECTION USING THE EXPRESSION WITH 1 / Vaa
        % if p.SDU == 1
        %     % compute second difference, Vaa
        %     Vdiff2 = zeros(nb, na, nz, ny);
            
        %     dF_adj = grids.a.dF;
        %     dF_adj(:,na) = grids.a.dB(:,na);
        %     grids_delta = dF_adj + grids.a.dB;
            
        %     V_i_plus_1_term = 2 * V(:,3:na,:,:) ./ (dF_adj(:,3:na) .* grids_delta(:,3:na));
        %     V_i_term = - 2 * V(:,2:na-1,:,:) .* (1./dF_adj(:,2:na-1) + 1./grids.a.dB(:,2:na-1))...
        %         ./ grids_delta(:,2:na-1);
        %     V_i_minus_1_term = 2 * V(:,1:na-2,:,:) ./ (grids.a.dB(:,1:na-2) .* grids_delta(:,1:na-2));
        %     Vdiff2(:,2:na-1,:,:) = V_i_plus_1_term + V_i_term + V_i_minus_1_term;
        %     Vdiff2(:,na,:,:) = (V(:,na-1,:,:) - V(:,na,:,:)) ./ (grids.a.dB(:,na) .^ 2);

        %     % compute adjustment term
        %     if p.invies == 1
        %         sdu_adj = 1 + (1-p.riskaver) * (Vdiff_SDU .^ 2) ./ Vdiff2;
        %     else
        %         sdu_adj = 1 + ((p.invies-p.riskaver) / (1-p.invies)) * (Vdiff_SDU .^2)...
        %             ./ (V .* Vdiff2);
        %     end
        %     sdu_adj(:,1,:,:) = 0;
        
        %     lowdiag = sdu_adj .* lowdiag;
        %     centdiag = sdu_adj .* centdiag;
        %     updiag = sdu_adj .* updiag;


        % end

        % lowdiag = circshift(reshape(lowdiag, nb*na, nz, ny), -nb);
        % updiag = circshift(reshape(updiag, nb*na, nz, ny), nb);

        % A = A + spdiags(centdiag(:),0,dim,dim);
        % A = A + spdiags(updiag(:),nb,dim,dim);
        % A = A + spdiags(lowdiag(:),-nb,dim,dim);

        A = A + compute_V_aa_terms(p, grids, nb, na, nz, ny);

        if p.SDU == 1

            VaB = zeros(nb, na, nz, ny);
            VaF = zeros(nb, na, nz, ny);
            
            VaB(:,2:na,:,:) = (V(:,2:na,:,:) - V(:,1:na-1,:,:)) ./ grids.a.dB(:,2:na);
            VaF(:,1:na-1,:,:) = (V(:,2:na,:,:) - V(:,1:na-1,:,:)) ./ grids.a.dF(:,1:na-1);
            Vdiff_SDU = (adriftB < 0) .* VaB + (adriftF > 0) .* VaF + ((adriftB>= 0) & (adriftF<=0)) .* Vdiff_SDU;
            % now add zeta(a) * Va^2 term
            if p.invies == 1
                adj_term = (1/2) * (grids.a.matrix * p.sigma_r) .^ 2 ...
                .* Vdiff_SDU * (1 - p.riskaver);
            else
                adj_term = (1/2) * (grids.a.matrix * p.sigma_r) .^ 2 ...
                    .* Vdiff_SDU ./ V * (p.invies - p.riskaver) / (1 - p.invies);
            end

            updiag = zeros(nb, na, nz, ny);
            lowdiag = zeros(nb, na, nz, ny);
            centdiag = zeros(nb, na, nz, ny);

            backward_delta = repmat(grids.a.dB, [1 1 nz ny]);
            forward_delta = repmat(grids.a.dF, [1 1 nz ny]);

            lowdiag(adriftB < 0) = - adj_term(adriftB < 0) ./ backward_delta(adriftB < 0);
            lowdiag(:,1,:,:) = 0;
            updiag(adriftF > 0) = adj_term(adriftF > 0) ./ forward_delta(adriftF > 0);
            updiag(:,na,:,:) = 0;
            centdiag = - lowdiag - updiag;

            stationary = (adriftB >= 0) & (adriftF <= 0);

            lowdiag = circshift(reshape(lowdiag, nb*na, nz, ny), -nb);
            updiag = circshift(reshape(updiag, nb*na, nz, ny), nb);

            A = A + spdiags(centdiag(:), 0, dim, dim);
            A = A + spdiags(updiag(:), nb, dim, dim);
            A = A + spdiags(lowdiag(:), -nb, dim, dim);
        end
    end

end

%% --------------------------------------------------------------------
% FUNCTIONS
% ---------------------------------------------------------------------
function Arisk_Vaa = compute_V_aa_terms(p, grids, nb, na, nz, ny)
    % This function computes the (1/2) * (a * sigma_r) ^2 * Vaa component
    % of the A matrix.

    top = false(nb, na);
    if p.OneAsset == 1
        asset_dB = repmat(grids.b.dB, [1 1 nz ny]);
        asset_dF = repmat(grids.b.dF, [1 1 nz ny]);

        top(nb, :, :, :) = true;

        risk_term = (grids.b.matrix * p.sigma_r) .^ 2;

        shift = 1;
    else
        asset_dB = repmat(grids.a.dB, [1 1 nz ny]);
        asset_dF = repmat(grids.a.dF, [1 1 nz ny]);

        top(:, na, :, :) = true;

        risk_term = (grids.a.matrix * p.sigma_r) .^ 2;

        shift = nb;
    end

    asset_dF(top) = asset_dB(top);
    asset_dSum = asset_dB + asset_dF;

    % (1/2) * (a * sigma_r) ^ 2 * V_{i+1} * (dF*(dB+dF)/2)
    V_i_plus_1_term = risk_term ./ (asset_dF .* asset_dSum);

    % - (1/2) * (a * sigma_r) ^ 2 * V_i * (1/dB + 1/dF) / ((dB+dF)/2)
    V_i_term = - risk_term .* (1./asset_dB + 1./asset_dF) ./ asset_dSum;
    
    % (1/2) * (a * sigma_r) ^ 2 * V_{i-1} * (dB*(dB+dF)/2)
    V_i_minus_1_term = risk_term ./ (asset_dB .* asset_dSum);

    % impose V' = 0 at the top of the grid
    V_i_term(top) = - risk_term(top) ./ (asset_dB(top) .* asset_dSum(top));
    V_i_plus_1_term(top) = 0;

    % shift entries to put on the A matrix diagonals in accordance with spdiags algorithm
    V_i_plus_1_term = circshift(reshape(V_i_plus_1_term, nb*na, nz, ny), shift);
    V_i_minus_1_term = circshift(reshape(V_i_minus_1_term, nb*na, nz, ny), -shift);
    dim = nb * na * nz * ny;

    Arisk_Vaa = spdiags(V_i_plus_1_term(:), shift, dim, dim)...
            + spdiags(V_i_term(:), 0, dim, dim)...
            + spdiags(V_i_minus_1_term(:), -shift, dim, dim);
end