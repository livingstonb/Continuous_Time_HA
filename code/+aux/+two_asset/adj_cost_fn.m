function adj_cost = adj_cost_fn(d,a_grid,p)
    % adjustment cost function chi(d)
    
    d_scaled = d./max(a_grid,p.a_lb);
    adj_cost = max(a_grid,p.a_lb) .* (p.chi0 * abs(d_scaled) + 1/(1+p.chi2) * (abs(d_scaled).^(1+p.chi2) * p.chi1^(-p.chi2)));

end