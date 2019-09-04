function P = penalty(a,penalty_param1,penalty_param2)
	P = zeros(size(a));
	P(a>=0) = 0;
	P(a<0) = penalty_param1 * abs(a(a<0)) .^ (penalty_param2);