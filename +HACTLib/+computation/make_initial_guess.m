function [V, gg] = make_initial_guess(p, grids, gridsKFE, income)
    % Makes initial guesses for the value function and distribution

    % Parameters
    % ----------
    % p : a Params object
    %
    % grids : a Grid object which holds the HJB grids
    %
    % gridsKFE : a Grid object which holds the KFE grids
    %
    % income : an Income object
    %
    % Returns
    % -------
    % V : value function guess, of shape (nb, na, nz, ny)
    %
    % gg : distribution guess, of shape (nb_KFE, na_KFE, nz, ny)

	nb = p.nb; na = p.na; 
	nz = p.nz; ny = income.ny;
	dim = nb * na * nz * ny;

	import HACTLib.aux.sparse_diags

    if numel(p.rhos) > 1
        rho_mat = reshape(p.rhos,[1 1 numel(p.rhos) 1]);
        rho_mat = repmat(rho_mat, [nb na 1 ny]);
        rho_mat = sparse_diags(rho_mat(:), 0);
    else
        rho_mat = p.rho * speye(nb*na*nz*ny);
    end
    rho_mat = rho_mat + p.deathrate * speye(nb*na*nz*ny);


    %% --------------------------------------------------------------------
    % GUESS FOR VALUE FUNCTION
    % ---------------------------------------------------------------------
	% liquid returns grid
	r_b_mat = p.r_b .* (grids.b.matrix>=0) +  p.r_b_borr .* (grids.b.matrix<0);

	% consumption guess
	r_b_mat_adj = r_b_mat;
	r_b_mat_adj(r_b_mat<=0.001) = 0.001; % mostly for r_b <= 0, enforce a reasonable guess
    r_a_adj = (p.r_a <= 0.001) * 0.001 + (p.r_a > 0.001) * p.r_a;
	c_0 = (1-p.directdeposit - p.wagetax) * income.y.matrix ...
            + (r_a_adj + p.deathrate*p.perfectannuities) * grids.a.matrix...
            + (r_b_mat_adj + p.deathrate*p.perfectannuities) .* grids.b.matrix + p.transfer;

    if p.SDU
        u = p.rho * HACTLib.aux.u_fn(c_0, p.invies);
    else
        u = HACTLib.aux.u_fn(c_0, p.riskaver);
    end

    inctrans = income.full_income_transition_matrix(p, u);
    
    if p.sigma_r > 0
        % Vaa term
        deltas = grids.a.dB + grids.a.dF;
        deltas(:, 1) = 2 * grids.a.dF(:, 1);
        deltas(:, na) = 2 * grids.a.dB(:, na);

        updiag = repmat(1 ./ grids.a.dF, [1 1 nz ny]);
        updiag(:,na,:,:) = repmat(1 ./ grids.a.dB(:,na), [1 1 nz ny]);

        centdiag = - repmat(1 ./ grids.a.dF + 1 ./ grids.a.dB, [1 1 nz ny]);
        centdiag(:,1,:,:) = repmat(1 ./ grids.a.dF(:,1), [1 1 nz ny]);
        centdiag(:,na,:,:) = -repmat(1 ./ grids.a.dB(:,na), [1 1 nz ny]);

        lowdiag = repmat(1 ./ grids.a.dB, [1 1 nz ny]);
        lowdiag(:,1,:,:) = -repmat(1 ./ grids.a.dF(:,1), [1 1 nz ny]);
        
        risk_adj = (grids.a.matrix .* p.sigma_r) .^ 2;

        updiag = risk_adj .* updiag ./ deltas;
        centdiag = risk_adj .* centdiag ./ deltas;
        lowdiag = risk_adj .* lowdiag ./ deltas;

        Arisk = HACTLib.aux.sparse_diags(...
            [lowdiag(:), centdiag(:), updiag(:)],...
            [-nb, 0, nb]);
    else
        Arisk = sparse(dim, dim);
    end

    V = (rho_mat - inctrans - Arisk) \ u(:);
    V = reshape(V, nb, na, nz, ny);

    %% --------------------------------------------------------------------
    % GUESS FOR EQUILIBRIUM DISTRIBUTION
    % ---------------------------------------------------------------------
    gg0 = ones(p.nb_KFE,p.na_KFE,p.nz,income.ny);
    gg0 = gg0 .* permute(repmat(income.ydist,[1 p.nb_KFE p.na_KFE p.nz]),[2 3 4 1]);
    if p.OneAsset
        gg0(:,gridsKFE.a.vec>0,:,:) = 0;
    end
    gg0 = gg0 / sum(gg0(:));
    gg0 = gg0 ./ gridsKFE.trapezoidal.matrix;
    gg = gg0;

end
