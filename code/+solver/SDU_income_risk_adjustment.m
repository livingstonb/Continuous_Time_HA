function ez_adj = SDU_income_risk_adjustment(p, Vn, income)
    % computes the risk adjustment in income transition rates
    % when households have stochastic differential utility

    assert(p.SDU == 1, "This function should not be called with SDU off.");

    
    shape = size(Vn);
    nb = shape(1);
    na = shape(2);
    
    nz = p.nz;
    ny = income.ny;
    
    if p.invies ~= 1
        ez_adj_0 = reshape(Vn, na*nb*nz, 1, ny) ./ reshape(Vn, na*nb*nz, ny, 1);
        ez_adj_1 = ((1-p.invies) ./ (1-p.riskaver))...
            .* ( (ez_adj_0 .^ ((1-p.riskaver)./(1-p.invies)) - 1) ./ (ez_adj_0 - 1) );
        
    else
        ez_adj_0 = (1-p.riskaver) * (reshape(Vn, na*nb*nz, 1, ny) ...
            - reshape(Vn, na*nb*nz, ny, 1));
        ez_adj_1 = (exp(ez_adj_0) - 1) ./ (ez_adj_0);
    end
    
    ez_adj = ez_adj_1 .* shiftdim(income.ytrans, -1);

    for kk = 1:income.ny
        idx_k = ~ismember(1:ny, kk);
        ez_adj(:,kk,kk) = -sum(ez_adj(:, kk, idx_k), 3);
    end
end