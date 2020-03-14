function z = interpolate_cdf(pmf, values, vrange,...
	query_point, porder)

	tmp = sortrows([values, pmf]);
	x = tmp(:,1);
	p = cumsum(tmp(:,2));

	dkeep = (x > vrange(1)) & (x < vrange(2));

	if any(dkeep)
		z = NaN;
	end

	x = x(dkeep);
	p = p(dkeep);

	poly_fitted = polyfit(x, p, porder);
	z = polyval(poly_fitted, query_point);
end