function Vn1 = solveHJB(p, A, income, Vn, u, nn, risk_adj)
    % This function updates the value function according to an
    % implicit or implicit-explicit updating scheme.

    % Parameters
    % ----------
    % p : a Params object which contains the parameters of the model
    %
    % income : an Income object which contains income variables
    %
    % A : the state transition matrix, sparse and of shape
    %   (nb*na*nz*ny, na*nb*nz*ny)
    %
    % Vn : the value function at this iteration, of shape
    %   (nb, na, nz, ny)
    %
    % u : the utility function at this iteration, of shape
    %   (nb, na, nz, ny)
    % nn : the given iteration number
    %
    % risk_adj : a risk adjustment array which comes from the
    %   Va^2 term in the SDU transformation with risky returns,
    %   of shape (nb, na, nz, ny)
    %
    % Returns
    % -------
    % Vn1 : the updated value function, of shape (nb, na, nz, ny)

    % grid lengths
    na = p.na;
    nb = p.nb;
    ny = numel(income.y.vec);
    nz = p.nz;

    if ~isempty(risk_adj)
        risk_adj_k = reshape(risk_adj,[],ny);
    end

    %% -----------------------------------------------------
    % RISK ADJUSTMENT FOR STOCHASTIC DIFFERENTIAL UTILITY
    % ------------------------------------------------------
    % the statement below returns [] for SDU = 0
    ez_adj = income.SDU_income_risk_adjustment(p, Vn);

    if p.implicit == 1
        %% -----------------------------------------------------
        % FULLY IMPLICIT UPDATING
        % ------------------------------------------------------
        % discount factor values
        if numel(p.rhos) > 1
            rhocol = repmat(kron(p.rhos', ones(nb*na,1)), ny, 1);
            rho_mat = spdiags(rhocol, 0, nb*na*nz*ny, nb*na*nz*ny);
        else
            rho_mat = p.rho * speye(nb*na*nz*ny);
        end

        assert(p.deathrate == 0, "Fully implicit assumes no death.")

        % add income transitions
        A = A + income.sparse_income_transitions(p, ez_adj, 'HJB');

        if isempty(risk_adj)
            RHS = p.delta_HJB * u(:) + Vn(:);
        else
            RHS = p.delta_HJB * u(:) + Vn(:) + p.delta_HJB * risk_adj(:);
        end
        
        % update V
        A = (rho_mat - A) * p.delta_HJB + speye(nb*na*nz*ny);
        Vn1 = A \ RHS;
        Vn1 = reshape(Vn1,nb,na,nz,ny);
    else
        %% -----------------------------------------------------
        % IMPLICIT-EXPLICIT UPDATING
        % ------------------------------------------------------
        
        u_k = reshape(u,[],ny);

        Vn1_k = NaN(nb*na*nz,ny);
        Vn_k = reshape(Vn,[],ny);
        Bik_all = cell(ny,1);

        % loop over income states
        for kk = 1:ny
            Bk = construct_Bk(p, income, kk, A, ez_adj);
            
            Bik_all{kk} = inverse(Bk);
            indx_k = ~ismember(1:ny,kk);

            if p.SDU == 1
                Vkp_stacked = sum(squeeze(ez_adj(:, kk, indx_k)) .* Vn_k(:, indx_k), 2);
            else
                Vkp_stacked = sum(repmat(income.ytrans(kk,indx_k),nb*na*nz,1) .* Vn_k(:,indx_k), 2);
            end

            qk = p.delta_HJB * u_k(:,kk) + Vn_k(:,kk) + p.delta_HJB*Vkp_stacked;
            if ~isempty(risk_adj)
                qk = qk + p.delta_HJB * risk_adj_k(:,kk);
            end
            
            Vn1_k(:,kk) = Bik_all{kk} * qk;
        end

        % Howard improvement step
        if (nn >= p.start_HIS) && (p.SDU == 0)
            Vn1_k = howard_improvement_step(p, income, ez_adj, Vn1_k, u_k, Bik_all);
        end
        Vn1 = reshape(Vn1_k,nb,na,nz,ny);
    end
end

function Bk = construct_Bk(p, income, k, A, ez_adj)
    % constructs the matrix Bk = (rho + deathrate - A)*delta + I
    % which serves as the divisor in the implicit-explicit update scheme

    % Parameters
    % ----------
    % p : a Params object
    %
    % income : an Income object
    %
    % k : index of the current income block within the A matrix
    %
    % A : the sparse transition matrix, of shape (nb*na*nz*ny, nb*na*nz*ny)
    %
    % ez_adj : the risk-adjusted income transitions, if utility is SDU, of
    %   shape (nb*na*nz, ny, ny)
    %
    % Returns
    % -------
    % Bk : the sparse matrix divisor in the implicit-explicit update scheme
    %   for income block k

    nb = p.nb;
    na = p.na;
    ny = income.ny;
    nz = p.nz;

    if numel(p.rhos) > 1
        rhocol = kron(p.rhos', ones(nb*na,1));
        rho_mat = spdiags(rhocol, nb*na*nz, nb*na*nz);
    else
        rho_mat = p.rho * speye(nb*na*nz);
    end

    indx_k = ~ismember(1:ny,k);
    Ak = A(1+(k-1)*(nb*na*nz):k*(nb*na*nz),1+(k-1)*(nb*na*nz):k*(nb*na*nz));
    Bk = p.delta_HJB * rho_mat + (1 + p.delta_HJB * p.deathrate) * speye(nb*na*nz)...
        - p.delta_HJB * Ak;

    if p.SDU == 1
        Bk = Bk - p.delta_HJB * spdiags(ez_adj(:, k, k), 0, nb*na*nz, nb*na*nz);
    else
        Bk = Bk - p.delta_HJB * income.ytrans(k, k) * speye(nb*na*nz);
    end

end

function Vn2_k = howard_improvement_step(p, income, ez_adj, Vn1_k, u_k, Bik_all)
    % technique to speed up convergence

    for jj = 1:p.maxit_HIS
        Vn2_k = NaN(p.nb*p.na*p.nz,income.ny);
        for kk = 1:income.ny
            indx_k = ~ismember(1:income.ny,kk);
            
            if p.SDU == 1
                Vkp_stacked = sum(squeeze(ez_adj(:, kk, indx_k)) .* Vn1_k(:, indx_k), 2);
            else
                Vkp_stacked = sum(repmat(income.ytrans(kk,indx_k),p.nb*p.na*p.nz,1) .* Vn1_k(:,indx_k),2);
            end
            qk = p.delta_HJB * u_k(:,kk) + Vn1_k(:,kk) + p.delta_HJB * Vkp_stacked;
            Vn2_k(:,kk) = Bik_all{kk} * qk;
        end

        dst = max(abs(Vn2_k(:) - Vn1_k(:)));
        Vn1_k = Vn2_k;
        if dst < p.crit_HIS
            break
        end
    end
end
