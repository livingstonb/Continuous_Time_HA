function d = opt_deposits(Va, Vb, a_grid, p)
	% computes the optimal deposit rate
	%
	% Parameters
	% ----------
	% Va : finite-difference approximation of the derivative of
	%	the value function with respect to illiquid assets, a
	%
	% Vb : finite-difference approximation of the derivative of
	%	the value function with respect to illiquid assets, b
	%
	% a_grid : illiquid asset levels
	%
	% p : a Params object
	%
	% Returns
	% -------
	% p : the optimal deposit rate over a_grid

	indx_0     = ((Va./Vb - 1 - p.chi0) <= 0) & ((Va./Vb - 1 + p.chi0) >= 0);
	indx_plus  = ((Va./Vb - 1 - p.chi0) > 0);
	indx_minus = ((Va./Vb - 1 + p.chi0) < 0);

	d = 0*indx_0 ...
	    + p.chi1 * (max(Va./Vb - 1 - p.chi0,0)).^(1/p.chi2) .* max(a_grid,p.a_lb) .* indx_plus ...
	    + (-p.chi1) * (max(-(Va./Vb - 1) - p.chi0,0)).^(1/p.chi2) .* max(a_grid,p.a_lb) .* indx_minus;

end