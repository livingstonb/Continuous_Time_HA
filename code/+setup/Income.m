classdef Income < handle
    % This class stores income grids and related variables, as well as useful
    % methods that act on income variables.
    
	properties (SetAccess = protected)
		% The income grid.
		y;

		% The log income grid.
		logy;
		
		% The square income transition matrix.
		ytrans;

		% The stationary distribution of the
		% income process.
		ydist;

		% The number of grid points for the
		% income process.
		ny;
	end

	methods
		function obj = Income(income_path, p, norisk)
  	 		% Class constructor, object can be initialized with an
  	 		% income path or passed no arguments if process is being
  	 		% set manually.
  	 		%
  	 		% Parameters
  	 		% ----------
  	 		% income_path : a string containing the path of the input
  	 		%	directory relative to the <root of this repository>/input
  	 		%
  	 		% p : An object containing at least the following fields:
  	 		%
  	 		%	nb, na, nb_KFE, na_KFE, and nz
  	 		%	- Grid sizes.
  	 		%
  	 		%	MPL
  	 		%	- Mean efficiency of labor. The income process will be
  	 		%	normalized so that mean income is MPL / 4.
  	 		%
  	 		% norisk : A flag which is true if income risk is turned off,
  	 		%	in which case income is equal to the mean for all households.
  	 		%	Default = false.

  	 		if nargin > 0
	  	 		if norisk
	  	 			obj.ydist = 1;
	  	 			obj.ytrans = 0;

	  	 			obj.y.vec = 1/4;
	  	 			obj.logy.vec = log(obj.y.vec);
	  	 		else
		            logy = load(fullfile(income_path, 'ygrid_combined.txt'));
		            y = exp(logy);
		            
		            obj.ydist = load(fullfile(income_path, 'ydist_combined.txt'));
		            obj.ytrans = load(fullfile(income_path, 'ymarkov_combined.txt'));
		            obj.ydist = aux.stat_dist(obj.ytrans');
		            
					% normalize
					obj.y.vec = y ./ (y' * obj.ydist * 4);
					obj.y.vec = p.MPL * obj.y.vec;
					
					obj.logy.vec = log(obj.y.vec);

					obj.generate(p);
		        end
	        end
        end

        function set_ydist(obj, ydist)
        	% Manually set the stationary distribution of income.
        	assert(isvector(ydist), "ydist must be a vector")
        	obj.ydist = ydist(:);

        	obj.ny = numel(ydist);
        end

        function set_ygrid(obj, ygrid, MPL)
        	% Manually set the values on the income grid.

        	assert(isvector(ydist), "ygrid must be a vector.")
        	assert(~isempty(obj.ydist), "Must set ydist before ygrid.")
        	assert(numel(obj.ydist) == numel(ygrid),...
        		"ygrid and ydist have inconsistent shapes.")

        	% normalize to mean MPL / 4
        	obj.y.vec = ygrid(:) ./ (ygrid(:)' * obj.ydist * 4);
			obj.y.vec = MPL * obj.y.vec;
        end

        function set_ytrans(obj, ytrans)
        	% Manually set the income transition matrix.
        	assert(ismatrix(ytrans), "ytrans must be a square matrix.")
			assert(size(ytrans, 1) == size(ytrans, 2),...
				"ytrans must be a square matrix.")
			assert(~isempty(obj.ydist), "Must set ydist before ytrans.")
			assert(obj.ny == size(ytrans, 1),...
				"ytrans and ydist have inconsistent shapes.")

			obj.ytrans = ytrans;
        end

        function generate(obj, p)
        	% Generates useful income variables.

        	obj.ny = numel(obj.y.vec);
            obj.y.wide = reshape(obj.y.vec,[1 1 1 obj.ny]);
		    obj.y.matrix = repmat(obj.y.wide,[p.nb p.na p.nz 1]);
            obj.y.matrixKFE = repmat(obj.y.wide,[p.nb_KFE p.na_KFE p.nz 1]);
        end

        function inctrans = sparse_income_transitions(obj, p, ez_adj, modeltype)
		    % Generates a sparse, square matrix of income transition rates 
		    % for solving the HJB or KFE.
		    %
		    % Parameters
		    % ----------
		    % p : An object containing at least the following fields:
  	 		%
  	 		%	nb, na, nb_KFE, na_KFE, and nz
  	 		%	- Grid sizes.
  	 		%
  	 		%	SDU
  	 		%	- Boolean indicator for stochastic differential utility.
		    %
		    % ez_adj : An array of shape (nb*na*nz, ny, ny) which contains the
		    %	income transitions adjusted for stochastic differential
		    %	utility. This argument is required if SDU is used, otherwise
		    %	it is ignored.
		    %
		    % modeltype : Grid type indicator, 'HJB' or 'KFE'.
		    %
		    % Returns
		    % -------
		    % inctrans : A sparse matrix of income transition rates, risk adjusted
		    %	if utility is SDU, and of shape (nb*na*nz*ny, nb*na*nz*ny).

		    if strcmp(modeltype, 'HJB')
			    if ~p.SDU
			        % return exogenous income transition rates
			        inctrans = kron(obj.ytrans, speye(p.nb*p.na*p.nz));
			    else
			        % adjust according to SDU transformation
			        ix = repmat((1:p.na*p.nb*p.nz*obj.ny)', obj.ny, 1);
			        iy = repmat((1:p.na*p.nb*p.nz)', obj.ny*obj.ny, 1);
			        iy = iy + kron((0:obj.ny-1)', p.nb*p.na*p.nz*ones(p.nb*p.na*p.nz*obj.ny,1));
			        inctrans = sparse(ix, iy, ez_adj(:));
			    end
			elseif strcmp(modeltype, 'KFE')
				% return exogenous income transition rates
			    inctrans = kron(obj.ytrans, speye(p.nb_KFE*p.na_KFE*p.nz));
			end
		end

		function ez_adj = SDU_income_risk_adjustment(obj, p, Vn)
		    % Computes the risk-adjusted income transition rates
		    % when households have stochastic differential utility.
		    % Returns [] when utility is not SDU.
		    %
		    % Parameters
		    % ----------
		    % p : An object with at least the following required fields:
		    %
		    %		SDU
		    %		- Boolean indicator for the use of stochastic differential
		    %	  	  utility.
		    %
		    %		riskaver
		    %		- Coefficient of risk aversion.
		    %	
		    %	and if SDU is true, then p must also have the field:
		    %
		    %		invies
		    %		- The inverse of the intertemporal elasticity of substitution.
		    %
		    % Vn : the value function, of shape (nb, na, nz, ny)
		    %
		    % Returns
		    % -------
		    % ez_adj : if utility is not SDU, returns []. otherwise, returns
		    %	risk-adjusted income transitions, of shape (nb*na*nz, ny, ny)

		    if ~p.SDU
		        ez_adj = [];
		        return;
		    end

		    nz = p.nz;
		    ny = obj.ny;
		    
		    shape = size(Vn);
		    nb = shape(1);
		    na = shape(2);
		    
		    if p.invies ~= 1
		    	if p.riskaver == 1
		    		error("Riskaver = 1, IES ~= 1 is not supported")
		    	end

		        ez_adj_0 = reshape(Vn, na*nb*nz, 1, ny) ./ reshape(Vn, na*nb*nz, ny, 1);
		        ez_adj_1 = ((1-p.invies) ./ (1-p.riskaver))...
		            .* ( (ez_adj_0 .^ ((1-p.riskaver)./(1-p.invies)) - 1) ./ (ez_adj_0 - 1) );
		    else
		        ez_adj_0 = (1-p.riskaver) * (reshape(Vn, na*nb*nz, 1, ny) ...
		            - reshape(Vn, na*nb*nz, ny, 1));
		        ez_adj_1 = (exp(ez_adj_0) - 1) ./ (ez_adj_0);
		    end
		    
		    ez_adj = ez_adj_1 .* shiftdim(obj.ytrans, -1);

		    for kk = 1:obj.ny
		        idx_k = ~ismember(1:ny, kk);
		        ez_adj(:,kk,kk) = -sum(ez_adj(:, kk, idx_k), 3);
		    end
		end
	end
end