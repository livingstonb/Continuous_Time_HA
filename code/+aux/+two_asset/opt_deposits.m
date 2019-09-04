function d = opt_deposits(Va,Vb,a_grid,p)

indx_0     = ((Va./Vb - 1 - p.chi0) <= 0) & ((Va./Vb - 1 + p.chi0) >= 0);
indx_plus  = ((Va./Vb - 1 - p.chi0) > 0);
indx_minus = ((Va./Vb - 1 + p.chi0) < 0);

d = 0*indx_0 ...
    + p.chi1 * (max(Va./Vb - 1 - p.chi0,0)).^(1/p.chi2) .* max(a_grid,p.a_lb) .* indx_plus ...
    + (-p.chi1) * (max(-(Va./Vb - 1) - p.chi0,0)).^(1/p.chi2) .* max(a_grid,p.a_lb) .* indx_minus;

end