function x = ra_calibrator(r_a, runopts, p)
	% This function solves the model for a given
	% value of illiquid returns

	% Set new illiquid returns
	p.reset_returns(p.r_b, exp(r_a));

    % Solve model
	stats = main_two_asset(runopts, p);

	% Compute distance from target
	x = stats.liqw - 0.5;

    fprintf("For r_a = %f:\n", p.r_a);
	fprintf("Liquid wealth = %f\n", stats.liqw);
end