function x = risk_premium_calibrator(returns, runopts, p)
	% This function solves the model for given values of
	% liquid returns and illiquid returns

	% Set new values for returns
	p.reset_returns(exp(returns(1))-0.05, exp(returns(2)));
	
	% Solve model
	stats = main_two_asset(runopts, p);

	% Compute distance from target
	x = stats.liqw - 0.5;
	x(2) = stats.totw - p.targetAY;

    fprintf("For r_b = %f", p.r_b);
	fprintf(" and r_a = %f:\n", p.r_a);
	fprintf("Liquid wealth = %f\n", stats.liqw);
	fprintf("Total wealth = %f\n", stats.totw);

end