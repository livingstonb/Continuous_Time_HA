classdef Income < handle
    % This class stores income grids and related variables, as well as useful
    % methods that act on income variables

    % The income process is read from files in the specified directory
    % and normalized so that mean quarterly income = 1/4
    
	properties (SetAccess = protected)
		logy;
		y;
		ytrans;
		ydist;
		ny;
	end

	methods
		function obj = Income(income_path, p, dimsHJB, dimsKFE, norisk)
  	 		% creates the income grids

  	 		% Parameters
  	 		% ----------
  	 		% income_path : a string containing the path of the input
  	 		%	directory
  	 		%
  	 		% p : a Params object which contains the parameters of the model
  	 		%
  	 		% dimsHJB : a row vector containing the sizes of the first four
  	 		%	dimensions for the HJB grids (e.g. [nb, na, nz, ny])
  	 		%
  	 		% dimsKFE : a row vector containing the sizes of the first four
  	 		%	dimensions for the KFE grids
  	 		%
  	 		% norisk : a flag which is true if income risk is turned off,
  	 		%	in which case income is equal to the mean for all households
  	 		%
  	 		% Returns
  	 		% -------
  	 		% an Income object

  	 		if norisk
  	 			obj.ydist = 1;
  	 			obj.ytrans = 0;

  	 			obj.y.vec = 1/4;
  	 			obj.logy.vec = log(obj.y.vec);
  	 		else
	            incdir = [income_path p.DirIncomeProcess '/'];
	            logy = load([incdir 'ygrid_combined.txt']);
	            y = exp(logy);
	            
	            obj.ydist = load([incdir 'ydist_combined.txt']);
	            obj.ytrans = load([incdir 'ymarkov_combined.txt']);
	            obj.ydist = aux.stat_dist(obj.ytrans');
	            
				% normalize
				obj.y.vec = y ./ (y' * obj.ydist * 4);
				obj.y.vec = p.MPL * obj.y.vec;
				
				obj.logy.vec = log(obj.y.vec);
	        end

	        obj.ny = numel(obj.y.vec);
            obj.y.wide = reshape(obj.y.vec,[1 1 1 obj.ny]);
		    obj.y.matrix = repmat(obj.y.wide,[dimsHJB 1]);
            obj.y.matrixKFE = repmat(obj.y.wide,[dimsKFE 1]);
        end

        function inctrans = sparse_income_transitions(obj, p, ez_adj)
		    % generate a sparse matrix of income transitions for solve 
		    % the HJB or KFE

		    % Parameters
		    % ----------
		    % p : a Params object which contains the parameters of the model
		    %
		    % ez_adj : an array of shape (nb*na*nz, ny, ny) which contains the
		    %	income transitions adjusted for stochastic diff utility
		    %
		    % Returns
		    % -------
		    % inctrans : a sparse matrix of income transition rates, risk adjusted
		    %	if utility is SDU, and of shape (nb*na*nz*ny, nb*na*nz*ny)

		    if p.SDU == 0
		        % return exogenous income transition rates
		        inctrans = kron(obj.ytrans, speye(p.nb*p.na*p.nz));
		    else
		        % adjust according to SDU transformation
		        ix = repmat((1:p.na*p.nb*p.nz*obj.ny)', obj.ny, 1);
		        iy = repmat((1:p.na*p.nb*p.nz)', obj.ny*obj.ny, 1);
		        iy = iy + kron((0:obj.ny-1)', p.nb*p.na*p.nz*ones(p.nb*p.na*p.nz*obj.ny,1));
		        inctrans = sparse(ix, iy, ez_adj(:));
		    end
		end

		function ez_adj = SDU_income_risk_adjustment(obj, p, Vn)
		    % computes the risk adjustment in income transition rates
		    % when households have stochastic differential utility

		    % Parameters
		    % ----------
		    % p : a Params object which contains the parameters of the model
		    %
		    % Vn : the value function, of shape (nb, na, nz, ny)
		    %
		    % Returns
		    % -------
		    % ez_adj : if utility is not SDU, returns []. otherwise, returns
		    %	risk-adjusted income transitions, of shape (nb*na*nz, ny, ny)

		    if p.SDU ~= 1
		        ez_adj = [];
		        return;
		    end
		    
		    shape = size(Vn);
		    nb = shape(1);
		    na = shape(2);
		    
		    nz = p.nz;
		    ny = obj.ny;
		    
		    if p.invies ~= 1
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