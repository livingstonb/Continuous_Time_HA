function fraction_below = find_constrained(values,pmf,thresholds)
	% given values and a pmf, this functions finds the fraction
	% of households below given thresholds
	sorted_by_values = sortrows([values(:) pmf(:)]);
	values_sorted = sorted_by_values(:,1);
	cdf_sorted = cumsum(sorted_by_values(:,2));

	pmf_support = sorted_by_values(:,2) > 1e-5;
	values_sorted = values_sorted(pmf_support);
	cdf_sorted = cdf_sorted(pmf_support);

	[values_unique,unique_indices] = unique(values_sorted,'last');
	cdf_unique = cdf_sorted(unique_indices);

	values_interp = griddedInterpolant(values_unique,cdf_unique,'linear');
	fraction_below = values_interp(thresholds);
end