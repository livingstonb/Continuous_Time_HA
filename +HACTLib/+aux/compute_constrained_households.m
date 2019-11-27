function fraction_lt = compute_constrained_households(values, pmf, thresholds)
	% Approximates the fraction of households for which
	% some variable is less than or equal to certain
	% thresholds, using interpolation.
	%
	% Parameters
	% ----------
	% values : Values over the state space.
	%
	% pmf : Probabilities associated with the
	%	entries in values.
	%
	% thresholds : A vector containing the
	%	the thresholds to compute.
	%
	% Returns
	% -------
	% fraction_lt : The fraction of households
	% 	less than or equal to each of the thresholds.

	validate_inputs(values, pmf, thresholds);

	sorted_by_values = sortrows([values(:) pmf(:)]);
	values_sorted = sorted_by_values(:,1);
	cdf_sorted = cumsum(sorted_by_values(:,2));

	[values_unique,unique_indices] = unique(values_sorted, 'last');
	cdf_unique = cdf_sorted(unique_indices);

	values_interp = griddedInterpolant(...
		values_unique,cdf_unique,' linear', 'linear');

	pmf_vec = pmf(:);
	fraction_lt = values_interp(thresholds(:));

	for ii = 1:numel(thresholds)
		if thresholds(ii) == min(values(:))
			% Threshold equal to minimum value, compute
			% probability of being in that state.
			fraction_lt(ii) = sum(pmf_vec(values(:)==min(values(:))));
		elseif thresholds(ii) < min(values(:))
			fraction_lt = 0;
		elseif thresholds(ii) >= max(values(:));
			fraction_lt(ii) = 1;
		end
	end
end

function validate_inputs(values, pmf, thresholds)
	assert(numel(values(:)) == numel(pmf(:)),...
		"Values and probabilities have inconsistent shapes");
	assert(~(isempty(values) || isempty(pmf) || isempty(thresholds)),...
		"One or more imputs is empty");
end