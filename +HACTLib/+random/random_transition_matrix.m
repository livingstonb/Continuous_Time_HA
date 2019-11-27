function A = random_transition_matrix(nb, na, varargin)
	% Constructs a random, sparse,
	% transition matrix with row sums of zero. .
	%
	% Inputs
	% ------
	% nb : Number of points on the liquid asset grid.
	%
	% na : Number of points on the illiquid asset grid.
	%	Must be >= 1.
	%
	% ny : Number of points on the income grid. Must
	%	be >= 1.
	%
	% Outputs
	% -------
	% A : A sparse transition matrix. 'A' will have
	%	a sparsity structure consistent with that
	%	found in an HACT model, given grid sizes of nb,
	%	na, and ny.

	options = parse_keyword_inputs(varargin{:});
	n_states = nb * na * options.ny;

	bdrift = (rand(nb, na, options.ny) - 0.5) * options.max_bdrift;
	bdrift_pos = max(bdrift, 0);
    bdrift_pos(nb, :, :) = 0;
    bdrift_pos = bdrift_pos(:);
    
	bdrift_neg = -min(bdrift, 0);
    bdrift_neg(1, :, :) = 0;
    bdrift_neg = bdrift_neg(:);
    bdrift_center = -bdrift_pos - bdrift_neg;
    
    bdrift_pos = [0; bdrift_pos(1:end-1)];
	bdrift_neg = [bdrift_neg(2:end); 0];

	A = spdiags(bdrift_pos, 1, n_states, n_states)...
		+ spdiags(bdrift_center, 0, n_states, n_states)...
		+ spdiags(bdrift_neg, -1, n_states, n_states);

	adrift = (rand(nb, na, options.ny) - 0.5) * options.max_adrift;
    adrift_pos = max(adrift, 0);
    adrift_pos(:,na,:) = 0;
    adrift_pos = adrift_pos(:);
    
	adrift_neg = -min(adrift, 0);
    adrift_neg(:,1,:) = 0;
    adrift_neg = adrift_neg(:);
    adrift_center = -adrift_pos - adrift_neg;
    
    adrift_pos = [zeros(nb, 1); adrift_pos(1:end-nb)];
	adrift_neg = [adrift_neg(nb+1:end); zeros(nb, 1)];

	A = A + spdiags(adrift_pos, nb, n_states, n_states)...
		+ spdiags(adrift_center, 0, n_states, n_states)...
		+ spdiags(adrift_neg, -nb, n_states, n_states);

	if nargin <= 2
		return
	end

	ytrans = rand(options.ny, options.ny);
	ytrans = ytrans ./ sum(ytrans, 2);
	for k = 1:options.ny
		indx_k = ~ismember(1:options.ny, k);
		ytrans(k, k) = -sum(ytrans(k, indx_k), 2);
	end

	A = A + kron(ytrans, speye(nb*na));
end

function options = parse_keyword_inputs(varargin)
	parser = inputParser;
	addParameter(parser, 'ny', 1);
	addParameter(parser, 'max_bdrift', 50);
	addParameter(parser, 'max_adrift', 50);
	parse(parser, varargin{:});

	options = parser.Results;
	assert(options.ny >= 1,...
		"ny must be greater than or equal to 1")
	assert(options.max_bdrift > 0,...
		"max_bdrift must be positive")
	assert(options.max_adrift > 0,...
		"max_adrift must be positive")
end