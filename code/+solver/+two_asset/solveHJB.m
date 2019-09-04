function Vn1 = solveHJB(p,A,income,Vn,u,nn)
    % Updates the value function

    na = p.na;
    nb = p.nb;
    ny = numel(income.y.vec);
    nz = p.nz;
    
    if numel(p.rhos) > 1
        % switch y and z order to match A matrix
        rhocol = kron(p.rhos',ones(nb*na,1));
        rho_mat = spdiags(rhocol,0,nb*na*nz,nb*na*nz);
    else
        rho_mat = p.rho * speye(nb*na*nz);
	end

    % Update value function

    Vn1 = NaN(nb,na,ny,nz);
    Bik_all = cell(ny,1);
    for kk = 1:ny
        Ak = A(1+(kk-1)*(nb*na*nz):kk*(nb*na*nz),1+(kk-1)*(nb*na*nz):kk*(nb*na*nz));

        Bk = p.delta_HJB * rho_mat + ...
            (1 + p.delta_HJB * p.deathrate...
            - p.delta_HJB*income.ytrans(kk,kk))*speye(nb*na*nz)...
            - p.delta_HJB*Ak;
        
        Bik_all{kk} 	= inverse(Bk);
        uk_stacked 		= reshape(u(:,:,kk,:),nb*na*nz,1);
        Vk_stacked 		= reshape(Vn(:,:,kk,:),nb*na*nz,1);
        indx_k 			= ~ismember(1:ny,kk);
        Vn_permuted     = permute(Vn,[1 2 4 3]);
        Vkp_stacked 	= sum(repmat(income.ytrans(kk,indx_k),nb*na*nz,1) ...
                            .* reshape(Vn_permuted(:,:,:,indx_k),nb*na*nz,ny-1),2);
        qk 				= p.delta_HJB*uk_stacked + Vk_stacked + p.delta_HJB*Vkp_stacked;
        Vn1k_stacked 	= Bik_all{kk}*qk;
        Vn1(:,:,kk,:) 	= reshape(Vn1k_stacked,nb,na,1,nz);
    end

    % Howard improvement step
    if nn >= p.start_HIS
        for jj = 1:p.maxit_HIS
            Vn2 = NaN(nb,na,ny,nz);
            for kk = 1:ny
                uk_stacked 		= reshape(u(:,:,kk,:),nb*na*nz,1);
                Vk_stacked 		= reshape(Vn1(:,:,kk,:),nb*na*nz,1);
                indx_k 			= ~ismember(1:ny,kk);
                Vn1_permuted    = permute(Vn1,[1 2 4 3]);
                Vkp_stacked 	= sum(repmat(income.ytrans(kk,indx_k),nb*na*nz,1) ...
                                    .* reshape(Vn1_permuted(:,:,:,indx_k),nb*na*nz,ny-1),2);
                qk 				= p.delta_HJB*uk_stacked + Vk_stacked + p.delta_HJB*Vkp_stacked;
                Vn2k_stacked 	= Bik_all{kk}*qk;
                Vn2(:,:,kk,:) 		= reshape(Vn2k_stacked,nb,na,1,nz);
            end
            dst = max(abs(Vn2(:) - Vn1(:)));
            Vn1 = Vn2;
            if dst < p.crit_HIS
                break
            end
        end
    end

    
end