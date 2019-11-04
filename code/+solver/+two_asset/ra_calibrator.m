function x = ra_calibrator(r_a, runopts, p)

	p.reset_returns(p.r_b, exp(r_a));

    % run model
	stats = main_two_asset(runopts, p);

	x = stats.liqw - 0.5;

    fprintf("For r_a = %f:\n", p.r_a);
	fprintf("Liquid wealth = %f\n", stats.liqw);
end