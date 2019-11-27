function gini = compute_gini(level, distr)
	% Computes the approximate Gini index.
	%
	% Parameters
	% ----------
	% level : Levels of the desired variable (e.g. wealth).
	%
	% distr : Probabilities associated with the entries
	%	in 'levels'.
	%
	% Returns
	% -------
	% gini : The Gini index, a scalar.

	validate_inputs(level, distr);

    % sort distribution and levels by levels
    sort1 = sortrows([level(:),distr(:)]);
    level_sort = sort1(:,1);
    dist_sort  = sort1(:,2);
    S = [0;cumsum(dist_sort .* level_sort)];
    gini = 1 - dist_sort' * (S(1:end-1)+S(2:end)) / S(end);
end

function validate_inputs(level, distr)
	assert(numel(level(:)) == numel(distr(:)),...
		"Inputs have inconsistent sizes");
	assert(~isempty(level),...
		"One or both inputs is empty");
end