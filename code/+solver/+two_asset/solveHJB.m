function Vn1 = solveHJB(p,A,income,Vn,u,nn)
    % Updates the value function

    na = p.na;
    nb = p.nb;
    ny = numel(income.y.vec);
    nz = p.nz;

    if (p.SDU == 1) && (p.invies ~= 1)
        ez_adj_0 = reshape(Vn, na*nb*nz, 1, ny) ./ reshape(Vn, na*nb*nz, ny, 1);
        ez_adj_1 = ((1-p.invies) ./ (1-p.riskaver))...
            .* ( (ez_adj_0 .^ ((1-p.riskaver)./(1-p.invies)) - 1) ./ (ez_adj_0 - 1) );
        
    elseif (p.SDU == 1) && (p.invies == 1)
        ez_adj_0 = (1-p.riskaver) * (reshape(Vn, na*nb*nz, 1, ny) - reshape(Vn, na*nb*nz, ny, 1));
        ez_adj_1 = (exp(ez_adj_0) - 1) ./ (ez_adj_0);
    end
    
    if p.SDU == 1
        ez_adj = ez_adj_1 .* shiftdim(income.ytrans, -1);

        for kk = 1:ny
            idx_k = ~ismember(1:ny, kk);
            ez_adj(:,kk,kk) = -sum(ez_adj(:, kk, idx_k), 3);
        end
    end

    if p.SDU == 0
        inctrans = kron(income.ytrans, speye(nb*na*nz));
    else
        ix = repmat((1:na*nb*nz*ny)', ny, 1);
        iy = repmat((1:na*nb*nz)', ny*ny, 1);
        iy = iy + kron((0:ny-1)', na*nb*nz*(ones(na*nb*nz*ny,1)));
        inctrans = sparse(ix, iy, ez_adj(:));
    end

    
   if (numel(p.rhos) > 1) && (p.implicit == 1)
        rhocol = repmat(kron(p.rhos', ones(nb*na,1)), ny, 1);
        rho_mat = spdiags(rhocol, 0, nb*na*nz*ny, nb*na*nz*ny);
    elseif p.implicit == 1
        rho_mat = p.rho * speye(nb*na*nz*ny);
    end

    %% -----------------------------------------------------
    % FULLY IMPLICIT UPDATING
    % ------------------------------------------------------
    if p.implicit == 1

        assert(p.deathrate == 0, "Fully implicit assumes no death.")

        % add income transitions
        if p.SDU == 0
            inctrans = kron(income.ytrans, speye(nb*na*nz));
        else
            ix = repmat((1:na*nb*nz*ny)', ny, 1);
            iy = repmat((1:na*nb*nz)', ny*ny, 1);
            iy = iy + kron((0:ny-1)', nb*na*nz*ones(nb*na*nz*ny,1));
            inctrans = sparse(ix, iy, ez_adj(:));
        end

        A = A + inctrans;
        clear inctrans ix iy ez_adj;
        RHS = p.delta_HJB * u(:) + Vn(:);
        
        A = (rho_mat - A) * p.delta_HJB + speye(nb*na*nz*ny);
        % [L,U] = ilu(B,struct('type','ilutp','droptol',1e-4));
        % Vn1 = ((repmat(rho_mat, ny, ny) - A) * p.delta_HJB + speye(nb*na*nz*ny)) \ RHS;
        % Vn1 = bicgstab(B, RHS, 1e-7, 10000, [], [], Vn(:));
        % Vn1 = gmres(B, RHS, 50, 1e-7, 300, [], [], Vn(:));
        Vn1 = A \ RHS;
        Vn1 = reshape(Vn1, nb, na, nz, ny);

    %% -----------------------------------------------------
    % IMPLICIT-EXPLICIT UPDATING
    % ------------------------------------------------------
    else
        % Update value function
        u_k = reshape(u,[],ny);

        Vn1_k = NaN(nb*na*nz,ny);
        Vn_k = reshape(Vn,[],ny);
        Bik_all = cell(ny,1);
        for kk = 1:ny
            Ak = A(1+(kk-1)*(nb*na*nz):kk*(nb*na*nz),1+(kk-1)*(nb*na*nz):kk*(nb*na*nz));

            Bk = p.delta_HJB * rho_mat + ...
                    (1 + p.delta_HJB * p.deathrate) * speye(nb*na*nz)...
                    - p.delta_HJB*Ak;
            if p.SDU == 1
                Bk = Bk - p.delta_HJB * spdiags(ez_adj(:, kk, kk), 0, nb*na*nz, nb*na*nz);
            else
                Bk = Bk - p.delta_HJB * income.ytrans(kk, kk) * speye(nb*na*nz);
            end
            
            Bik_all{kk} 	= inverse(Bk);
            indx_k 			= ~ismember(1:ny,kk);

            if p.SDU == 1
                Vkp_stacked = sum(squeeze(ez_adj(:, kk, indx_k)) .* Vn_k(:, indx_k), 2);
            else
                Vkp_stacked = sum(repmat(income.ytrans(kk,indx_k),nb*na*nz,1) .* Vn_k(:,indx_k), 2);
            end

            qk = p.delta_HJB * u_k(:,kk) + Vn_k(:,kk) + p.delta_HJB*Vkp_stacked;
            Vn1_k(:,kk) = Bik_all{kk} * qk;
        end

        % Howard improvement step
        if (nn >= p.start_HIS) && (p.SDU==0)
            for jj = 1:p.maxit_HIS
                Vn2_k = NaN(nb*na*nz,ny);
                for kk = 1:ny
                    indx_k 			= ~ismember(1:ny,kk);
                    Vkp_stacked 	= sum(repmat(income.ytrans(kk,indx_k),nb*na*nz,1) .* Vn1_k(:,indx_k),2);
                    qk 				= p.delta_HJB * u_k(:,kk) + Vn1_k(:,kk) + p.delta_HJB * Vkp_stacked;
                    Vn2_k(:,kk) 	= Bik_all{kk} * qk;
                end

                dst = max(abs(Vn2_k(:) - Vn1_k(:)));
                Vn1_k = Vn2_k;
                if dst < p.crit_HIS
                    break
                end
            end
        end

        Vn1 = reshape(Vn1_k,nb,na,nz,ny);
    end
end