function percentiles_values = compute_pct(values, pmf, percentiles)
	% Finds percentiles given a pmf over values.
	%
	% Parameters
	% ----------
	% values : Values at different states.
	%
	% pmf : Probability mass function over states,
	%	must have the same number of elements as
	%	'values'.
	%
	% percentiles : Percentiles requested, between 0 and 1,
	%	a row or column vector. May be a scalar.
	%
	% Returns
	% -------
	% percentiles_values : The values of the requested
	% 	percentiles, of the same shape as 'percentiles'.

	validate_inputs(values, pmf, percentiles);
	
	sorted_by_values = sortrows([values(:) pmf(:)]);
	values_sorted = sorted_by_values(:,1);
	cdf_sorted = cumsum(sorted_by_values(:,2));

	[cdf_unique,unique_indices] = unique(cdf_sorted,'last');
	values_sorted_unique = values_sorted(unique_indices);

	cdf_interp = griddedInterpolant(cdf_unique, values_sorted_unique, 'linear');
	percentiles_values = cdf_interp(percentiles(:));
	percentiles_values = reshape(percentiles_values, size(percentiles));
end

function validate_inputs(values, pmf, percentiles)
	assert(numel(values(:)) == numel(pmf(:)),...
		"Values and probabilities have inconsistent shapes");
	assert(~(isempty(values) || isempty(pmf) || isempty(percentiles)),...
		"One or more imputs is empty");
	assert(all(percentiles(:) > 0) && all(percentiles(:) < 1),...
		"One of the requested percentiles is outside of (0,1)");
end