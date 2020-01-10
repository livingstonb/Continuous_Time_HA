function inctrans = full_income_transition_matrix(ytrans, dims, SDU, p, V)
    % Generates a sparse, square matrix of income transition rates 
    % over the entire state space. Applies a risk-adjustment when
    % utility is SDU.
    %
    % Required Inputs
    % ---------------
    % ytrans : The square income transition matrix. Rows
    %	should sum to 0.
    %
    % dims : A row vector of length three which contains the
    %	number of states along the non-income dimensions.
    %
    % Additional Required Inputs for SDU
    % ----------------------------------
    %
    % SDU : A boolean flag indicating whether or not utility
    %	is SDU. If SDU is false, then 'p' and 'V' are ignored.
    %
    % p : An object containing the attributes 'riskaver' and
    %	'invies'.
    %
    % V : The value function over the full state space.
    %
    % Returns
    % -------
    % inctrans : A sparse matrix of income transition rates,
    %	of shape (nb*na*nz*ny, nb*na*nz*ny). If utility is
    %	SDU, these transition rates are risk-adjusted.

    import HACTLib.Checks

    narginchk(2, 5);
    Checks.has_shape('income_transition_matrix', dims, [1, 3]);
    Checks.is_integer('income_transition_matrix', dims);
    nb = dims(1);
    na = dims(2);
    nz = dims(3);
    ny = size(ytrans, 1);

    if nargin == 2
    	inctrans = kron(ytrans, speye(nb*na*nz));
    elseif ~SDU
    	inctrans = kron(ytrans, speye(nb*na*nz));
    else
    	ez_adj = risk_adjusted_income_transitions_SDU(ytrans, dims, p, V);
	    ix = repmat((1:na*nb*nz*ny)', ny, 1);
	    iy = repmat((1:na*nb*nz)', ny*ny, 1);
	    iy = iy + kron((0:ny-1)', nb*na*nz*ones(nb*na*nz*ny, 1));
	    inctrans = sparse(ix, iy, ez_adj(:));
	end
end