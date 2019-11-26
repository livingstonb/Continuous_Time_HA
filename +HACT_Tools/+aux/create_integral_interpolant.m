function interpolant = create_integral_interpolant(grid_values, integrand_values, pmf)
	% Creates an interpolant that approximates the value of the integral
	% int_0^{epsilon} values(a)g(a)da for a given epsilon.
	%
	% Parameters
	% ----------
	% grid_values : Values at which the integrand is evaluated.
	%
	% integrand_values : Values of the integrand.
	%
	% pmf : The probability mass function over states.
	%
	% Results
	% -------
	% interpolant : A griddedInterpolant object such that interpolant(x)
	%	is the approximated value of the integral from 0 to x.

	validate_inputs(grid_values, integrand_values, pmf);

	sorted_inputs = sortrows([grid_values(:) integrand_values(:) pmf(:)]);
	grid_sorted = sorted_inputs(:,1);
	integrand_sorted = sorted_inputs(:,2);
	pmf_sorted = sorted_inputs(:,3);

	integral_values = cumsum(integrand_sorted .* pmf_sorted);

	[grid_unique, unique_inds] = unique(grid_sorted,'last');
	integral_unique = integral_values(unique_inds);

	interpolant = griddedInterpolant(...
		grid_unique, integral_unique, 'linear');
end

function validate_inputs(grid_values, integrand_values, pmf)
	assert(numel(grid_values(:)) == numel(integrand_values(:)),...
		"Inputs have inconsistent shapes");
	assert(numel(grid_values(:)) == numel(pmf(:)),...
		"Inputs have inconsistent shapes");
	assert(~(isempty(grid_values) || isempty(integrand_values) || isempty(pmf)),...
		"One or more imputs is empty");
end