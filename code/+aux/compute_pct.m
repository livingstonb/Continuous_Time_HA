function percentiles_values = compute_pct(values, pmf, percentiles)
	% finds percentiles, given a pmf over values

	% Parameters
	% ----------
	% values : values at states
	%
	% pmf : probability mass function over states
	%
	% percentiles : percentiles requested, between 0 and 1
	%
	% Returns
	% -------
	% percentiles_values : the values of the requested percentiles
	
	sorted_by_values = sortrows([values(:) pmf(:)]);
	values_sorted = sorted_by_values(:,1);
	cdf_sorted = cumsum(sorted_by_values(:,2));

	[cdf_unique,unique_indices] = unique(cdf_sorted,'last');
	values_sorted_unique = values_sorted(unique_indices);

	cdf_interp = griddedInterpolant(cdf_unique,values_sorted_unique,'linear');
	percentiles_values = cdf_interp(percentiles);
end