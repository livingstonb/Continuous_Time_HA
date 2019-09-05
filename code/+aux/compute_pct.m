function percentiles_values = compute_pct(values,pmf,percentiles)
	% finds percentiles, given a pmf over values

	% 'percentiles' input vector should be fractions, i.e.
	% to compute 99th percentile, use 0.99
	
	sorted_by_values = sortrows([values(:) pmf(:)]);
	values_sorted = sorted_by_values(:,1);
	cdf_sorted = cumsum(sorted_by_values(:,2));

	[cdf_unique,unique_indices] = unique(cdf_sorted,'last');
	values_sorted_unique = values_sorted(unique_indices);

	cdf_interp = griddedInterpolant(cdf_unique,values_sorted_unique,'linear');
	percentiles_values = cdf_interp(percentiles);
end