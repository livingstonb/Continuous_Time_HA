function inctrans = get_income_transitions(p, income, ez_adj)
    if p.SDU == 0
        % return exogenous income transition rates
        inctrans = kron(income.ytrans, speye(p.nb*p.na*p.nz));
    else
        % adjust according to SDU transformation
        ix = repmat((1:p.na*p.nb*p.nz*income.ny)', income.ny, 1);
        iy = repmat((1:p.na*p.nb*p.nz)', income.ny*income.ny, 1);
        iy = iy + kron((0:income.ny-1)', p.nb*p.na*p.nz*ones(p.nb*p.na*p.nz*income.ny,1));
        inctrans = sparse(ix, iy, ez_adj(:));
    end
end