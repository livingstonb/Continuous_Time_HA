function Vn1 = solveHJB(p,A,income,Vn,u,nn)
    % Updates the value function

    na = p.na;
    nb = p.nb;
    ny = numel(income.y.vec);
    nz = p.nz;
    
    if numel(p.rhos) > 1
        rhocol = kron(p.rhos',ones(nb*na,1));
        rho_mat = spdiags(rhocol,0,nb*na*nz,nb*na*nz);
    else
        rho_mat = p.rho * speye(nb*na*nz);
	end

    % Update value function
    u_k = reshape(u,[],ny);

    Vn1_k = NaN(nb*na*nz,ny);
    Vn_k = reshape(Vn,[],ny);
    Bik_all = cell(ny,1);
    for kk = 1:ny
        Ak = A(1+(kk-1)*(nb*na*nz):kk*(nb*na*nz),1+(kk-1)*(nb*na*nz):kk*(nb*na*nz));

        Bk = p.delta_HJB * rho_mat + ...
            (1 + p.delta_HJB * p.deathrate...
            - p.delta_HJB*income.ytrans(kk,kk)) * speye(nb*na*nz)...
            - p.delta_HJB*Ak;
        
        Bik_all{kk} 	= inverse(Bk);
        indx_k 			= ~ismember(1:ny,kk);
        Vkp_stacked 	= sum(repmat(income.ytrans(kk,indx_k),nb*na*nz,1) .* Vn_k(:,indx_k),2);
        qk 				= p.delta_HJB * u_k(:,kk) + Vn_k(:,kk) + p.delta_HJB*Vkp_stacked;
        Vn1_k(:,kk) 	= Bik_all{kk} * qk;
    end

    % Howard improvement step
    if nn >= p.start_HIS
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