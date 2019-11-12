function V = value_guess_risky_returns(p, grids, income)

	nb = p.nb;
	na = p.na;
	nz = p.nz;
	ny = income.ny;

	dim = nb * na * nz * ny;

	rho_mat = p.rho * speye(dim);

	% liquid returns grid
	r_b_mat = p.r_b .* (grids.b.matrix>=0) +  p.r_b_borr .* (grids.b.matrix<0);

	% consumption guess
	r_b_mat_adj = r_b_mat;
	r_b_mat_adj(r_b_mat<=0.001) = 0.001; % mostly for r_b <= 0, enforce a reasonable guess
	c_0 = (1-p.directdeposit - p.wagetax) * income.y.matrix ...
            + (p.r_a + p.deathrate*p.perfectannuities) * grids.a.matrix...
            + (r_b_mat_adj + p.deathrate*p.perfectannuities) .* grids.b.matrix + p.transfer;

    u = aux.u_fn(c_0, p.invies);

    if p.SDU == 1
        % risk-adjusted income transitions
        if p.invies ~= 1
            ez_adj_0 = reshape(u, na*nb*nz, 1, ny) ./ reshape(u, na*nb*nz, ny, 1);
            ez_adj_1 = ((1-p.invies) ./ (1-p.riskaver))...
                .* ( (ez_adj_0 .^ ((1-p.riskaver)./(1-p.invies)) - 1) ./ (ez_adj_0 - 1) );
            
        else
            ez_adj_0 = (1-p.riskaver) * (reshape(u, na*nb*nz, 1, ny) - reshape(u, na*nb*nz, ny, 1));
            ez_adj_1 = (exp(ez_adj_0) - 1) ./ (ez_adj_0);
        end
        
        ez_adj = ez_adj_1 .* shiftdim(income.ytrans, -1);

        for kk = 1:ny
            idx_k = ~ismember(1:ny, kk);
            ez_adj(:,kk,kk) = -sum(ez_adj(:, kk, idx_k), 3);
        end
        

        ix = repmat((1:na*nb*nz*ny)', ny, 1);
        iy = repmat((1:na*nb*nz)', ny*ny, 1);
        iy = iy + kron((0:ny-1)', nb*na*nz*ones(nb*na*nz*ny,1));
        inctrans = sparse(ix, iy, ez_adj(:));
    else
        inctrans = kron(income.ytrans, speye(nb * na * nz));
    end

    if p.sigma_r > 0
        % Vaa coeffs
        deltas = grids.a.dB + grids.a.dF;
        deltas(:, 1) = 2 * grids.a.dF(:, 1);
        deltas(:, na) = 2 * grids.a.dB(:, na);

        updiag = zeros(nb, na, nz, ny);
        centdiag = zeros(nb, na, nz, ny);
        lowdiag = zeros(nb, na, nz, ny);
        
        updiag(:, 1:na-1, :, :) = repmat(1 ./ grids.a.dF(:, 1:na-1), [1 1 nz ny]);

        centdiag(:, 1:na-1, :, :) = - repmat(1 ./ grids.a.dF(:, 1:na-1) + 1 ./ grids.a.dB(:, 1:na-1), [1 1 nz ny]);
        centdiag(:, na, :, :) = - repmat(1 ./ grids.a.dB(:, na), [1 1 nz ny]);

        lowdiag(:, 1:na-1, :, :) = repmat(1 ./ grids.a.dB(:, 1:na-1), [1 1 nz ny]);
        lowdiag(:, na, :, :) = repmat(1 ./ grids.a.dB(:, na), [1 1 nz ny]);
        
        updiag(:, 1, :, :) = 0;
        centdiag(:, 1, :, :) = 0;
        lowdiag(:, 1, :, :) = 0;

        risk_adj = (grids.a.matrix .* p.sigma_r) .^ 2;

        updiag = risk_adj .* updiag ./ deltas;
        centdiag = risk_adj .* centdiag ./ deltas;
        lowdiag = risk_adj .* lowdiag ./ deltas;

        updiag = circshift(reshape(updiag, nb*na, nz, ny), nb);
        lowdiag = circshift(reshape(lowdiag, nb*na, nz, ny), -nb);

        updiag = updiag(:);
        centdiag = centdiag(:);
        lowdiag = lowdiag(:);

        Arisk = spdiags(updiag, nb, dim, dim)...
        	+ spdiags(centdiag, 0, dim, dim)...
        	+ spdiags(lowdiag, -nb, dim, dim);
    else
        Arisk = sparse(dim, dim);
    end

    if p.SDU == 1
        V = (rho_mat - inctrans - Arisk) \ reshape(p.rho .* u, [], 1);
    else
        V = (rho_mat - inctrans - Arisk) \ u(:);
    end

end