function x = risk_premium_calibrator(returns, runopts, p)

	p.reset_returns(exp(returns(1))-0.01, exp(returns(2)));

	fprintf("Trying r_b = %f\n", p.r_b);
	fprintf("Trying r_a = %f\n", p.r_a);

	stats = main_two_asset(runopts, p);

	x = stats.liqw - 0.5;
	x(2) = stats.totw - 3.5;

	fprintf("Liquid wealth = %f\n", stats.liqw);
	fprintf("Total wealth = %f\n", stats.totw);

end