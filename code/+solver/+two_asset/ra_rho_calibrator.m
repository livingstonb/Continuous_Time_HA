function y = ra_rho_calibrator(x, runopts, p)
	% This function solves the model for a given
	% value of illiquid returns
    
    new_rho = 0.05 * abs(x(1)) / (1 + abs(x(1)));
    
    new_ra = p.r_b + 0.04 * abs(x(2)) / (1 + abs(x(2)));
    
    % Set new discount rate
    p.set("rho", new_rho);

	% Set new illiquid returns
	p.set("r_a", new_ra);

    % Solve model
	stats = main_two_asset(runopts, p);

	% Compute distance from target
    y = zeros(2, 1);
	y(1) = (stats.liqw - 0.5) ^ 2;
    y(2) = (stats.totw - p.targetAY) ^ 2;

    fprintf("For rho=%f:\n", p.rho);
    fprintf("For r_a = %f:\n", p.r_a);
	fprintf("Liquid wealth = %f\n", stats.liqw);
    fprintf("Total wealth = %f\n", stats.totw);
end