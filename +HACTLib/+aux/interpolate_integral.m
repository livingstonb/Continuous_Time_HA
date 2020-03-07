function interpolant = interpolate_integral(grid_values, integrand_values, pmf, varargin)
	% Creates an interpolant that approximates the value of the integral
	% int_0^{epsilon} values(a)g(a)da for a given epsilon.
	%
	% Parameters
	% ----------
	% gridValues : Values at which the integrand is evaluated.
	%
	% integrand_values : Values of the integrand.
	%
	% pmf : The probability mass function over states.
	%
	% Results
	% -------
	% interpolant : A griddedInterpolant object such that interpolant(x)
	%	is the approximated value of the integral from 0 to x.

	is_sorted = validate_inputs(...
		grid_values, integrand_values, pmf, varargin{:});

	if ~is_sorted
		sorted_inputs = sortrows([grid_values(:) integrand_values(:) pmf(:)]);
		grid_values = sorted_inputs(:,1);
		integrand_values = sorted_inputs(:,2);
		pmf = sorted_inputs(:,3);
	end

	integral_values = cumsum(integrand_values .* pmf);
	[grid_unique, unique_inds] = unique(grid_values, 'last');
	integral_unique = integral_values(unique_inds);

	interpolant = griddedInterpolant(...
		grid_unique, integral_unique, 'pchip', 'nearest');
end

function is_sorted = validate_inputs(grid_values, integrand_values, pmf , varargin)
	assert(numel(grid_values(:)) == numel(integrand_values(:)),...
		"Inputs have inconsistent shapes");
	assert(numel(grid_values(:)) == numel(pmf(:)),...
		"Inputs have inconsistent shapes");
	assert(~(isempty(grid_values) || isempty(integrand_values) || isempty(pmf)),...
		"One or more imputs is empty");

	is_sorted = false;
	if ~isempty(varargin)
		if islogical(varargin{1})
			is_sorted = varargin{1};
		end
	end
end