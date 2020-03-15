function r_b_out = net_liquid_returns(bgrid, r_b, r_b_borr)
	r_b_out = r_b .* (bgrid >=0 ) +  r_b_borr .* (bgrid < 0);
end