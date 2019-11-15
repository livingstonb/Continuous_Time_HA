function d = adj_cost_deriv_inverse(y,a_grid,p)
    % Inverse of the derivative of the adjustment cost function,
    % i.e. (chi'(d))^{-1}
	d = sign(y) .* max(a_grid,p.a_lb) .* p.chi1 .* (abs(y) - p.chi0).^(1/p.chi2);
end