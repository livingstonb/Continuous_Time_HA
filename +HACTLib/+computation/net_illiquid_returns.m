function r_a_net = net_illiquid_returns(avalues, amax, r_a, coeff1, coeff2)
	if nargin == 3
		coeff1 = 0.2;% 0.075;
		coeff2 = 5;% 9.0;
	end

	% r_a_net = r_a - 0.9 * r_a .* (avalues >= 0.4 * amax) .* (avalues ./ amax) .^ 2;
	r_a_net = r_a;

	% tax = coeff1 * (avalues ./ amax) .^ coeff2;
	% r_a_net =  r_a .* (1 - tax);
end