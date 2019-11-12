function chi_prime = adj_cost_deriv(d, a_grid, p)
	d_scaled = d ./ max(a_grid, p.a_lb);
	chi_prime = p.chi0 + sign(d) .* p.chi1 ^(-p.chi2) .* abs(d_scaled) .^ p.chi2;
end