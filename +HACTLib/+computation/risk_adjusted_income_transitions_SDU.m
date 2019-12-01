function ez_adj = risk_adjusted_income_transitions_SDU(ytrans, dims, p, V)
    % Computes the risk-adjusted income transition rates
    % when households have stochastic differential utility.
    % Returns [] when utility is not SDU.
    %
    % Parameters
    % ----------
    % ytrans : The square income transition matrix. Rows should
    %	sum to zero.
    %
    % dims : A row vector of length three containing the lengths
    %	of the non-income dimensions.
    %
    % p : An object with at least the following required fields:
    %
    %		riskaver
    %		- The coefficient of risk aversion.
    %
    %		invies
    %		- The inverse of the intertemporal elasticity of substitution.
    %
    % V : the value function, of shape (nb, na, nz, ny)
    %
    % Returns
    % -------
    % ez_adj : if utility is not SDU, returns []. otherwise, returns
    %	risk-adjusted income transitions, of shape (nb*na*nz, ny, ny)

    nb = dims(1);
    na = dims(2);
    nz = dims(3);
    ny = size(ytrans, 1);
    
    if p.invies ~= 1
    	if p.riskaver == 1
    		error("Riskaver = 1, IES ~= 1 is not supported")
    	end

        ez_adj_0 = reshape(V, na*nb*nz, 1, ny) ./ reshape(V, na*nb*nz, ny, 1);
        ez_adj_1 = ((1-p.invies) ./ (1-p.riskaver))...
            .* ( (ez_adj_0 .^ ((1-p.riskaver)./(1-p.invies)) - 1) ./ (ez_adj_0 - 1) );
    else
        ez_adj_0 = (1-p.riskaver) * (reshape(V, na*nb*nz, 1, ny) ...
            - reshape(V, na*nb*nz, ny, 1));
        ez_adj_1 = (exp(ez_adj_0) - 1) ./ (ez_adj_0);
    end
    
    ez_adj = ez_adj_1 .* shiftdim(ytrans, -1);

    for kk = 1:ny
        idx_k = ~ismember(1:ny, kk);
        ez_adj(:,kk,kk) = -sum(ez_adj(:, kk, idx_k), 3);
    end
end