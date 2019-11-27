classdef TransitionMatrixConstructor < handle
    % This class constructs the A matrix for the HJB or KFE
    % First, instantiate the class with the required arguments,
    % then call the construct method with the policy functions
    % and value function as arguments.

    properties (SetAccess = protected)
    	%	An object with the following required attributes:
    	%
    	%		deathrate
    	%		- The Poisson rate at which households die.
    	%
    	%		perfectannuities
    	%		- Boolean, true turns on annuities, false
    	%		  leaves them off.
    	%
    	%		directdeposit
    	%		- Fraction of income being paid out to the
    	%		  illiquid asset.
    	%
    	%		SDU
    	%		- A boolean indicator of whether or not stochastic
    	%		  differential utility is present.
    	%		
    	%		sigma_r
    	%		- The standard deviation of returns risk, typically
    	%		  zero.
    	%
    	%    and if SDU evaluates to true, then 'p' must also include:
    	%		
    	%		riskaver
    	%		- The coefficient of risk aversion.
    	%
    	%		invies
    	%		- The inverse of the intertemporal elasticity of
    	%		  substitution.
    	p;

    	% A Grid object.
        grids;

        nb;
        na;
        ny;
        nz;
        n_states;

        % Shape of the full state space, a row vector.
        shape;

        % returns risk is on or off
        returns_risk;

        % Grid deltas
        asset_dB;
        asset_dF;
        asset_dSum; % (deltaB + deltaF)

        % Boolean mask for top of the b (1-asset) or a (2-asset) grid
        top;

        % Pre-computed term when returns are risky
        risk_term;

        % Dimension of terms for upper/lower diagonals. Offset
        % from diagonal by 1 if shift_dimension = 1 and by 
        % nb if shift_dimension = 2.
        shift_dimension;

        % Matrix of income values.
        y_mat;

        % 'HJB' or 'KFE' grid
        gridtype;
    end

    methods
        function obj = TransitionMatrixConstructor(p, income, grids, returns_risk)
            % Class constructor.
            %
            % Parameters
            % ----------
            % p : An object with the required attributes outlined in the
            %	  Properties block.
            %
            % income : An object with the required attributes outlined in
            %	  the Properties block.
            %
            % grids : A Grid object.
            %
            % returns_risk : A boolean indicator for whether or not returns
            %	  risk is to be included in the transition matrix.

            obj.nb = numel(grids.b.vec);
            obj.na = numel(grids.a.vec);
            obj.nz = grids.nz;
            obj.ny = numel(income.y.vec);
            obj.n_states = obj.nb * obj.na * obj.nz * obj.ny;
            obj.shape = [obj.nb obj.na obj.nz obj.ny];

            obj.gridtype = grids.gtype;
            obj.grids = grids;
            obj.p = p;

            obj.returns_risk = returns_risk;

            if strcmp(obj.gridtype, 'KFE')
                obj.y_mat = income.y.matrixKFE;
            elseif strcmp(obj.gridtype, 'HJB')
                obj.y_mat = income.y.matrix;
            end

            if p.sigma_r > 0
                obj.perform_returns_risk_computations();
            end
        end

        %% --------------------------------------------------------------------
        % Construct the A Matrix
        % ---------------------------------------------------------------------
        function [A, stationary, drifts] = construct(obj, model, V)
            % Constructs the transition matrix.
            %
            % Parameters
            % ----------
            % model : A structure containing the policy functions, saving (s)
            %   and deposits (d) are used here.
            %
            % V : The value function, of shape (nb, na, nz, ny).
            %
            % Returns
            % -------
            % A : The sparse transition matrix of shape (nb*na*nz*ny, nb*na*nz*ny).
            %
            % stationary : A boolean mask indicating states where drift in the
            %   risky asset was neither backward nor forward. Risk computations
            %   for these states will be included outside the A matrix.
            
            drifts = obj.compute_asset_drifts(model);

            A = obj.compute_liquid_transitions(drifts.b_B, drifts.b_F);
            A = A + obj.compute_illiquid_transitions(drifts.a_B, drifts.a_F);

            stationary = [];
            if obj.returns_risk
            	A = A + obj.compute_Vaa_terms();
            end  
        endd
    end

    methods (Access=protected)
    	%% --------------------------------------------------------------------
        % Compute Asset Drifts
        % ---------------------------------------------------------------------
        function drifts = compute_asset_drifts(obj, model)
        	% Approximates the drifts of the liquid and illiquid assets.
        	% The input variable 'model' must contain the policy
        	% functions 's' and 'd'.

            if strcmp(obj.gridtype,'KFE')
                drifts.a_B = min(model.d + obj.grids.a.matrix * (obj.p.r_a + obj.p.deathrate*obj.p.perfectannuities)...
                            + obj.p.directdeposit * obj.y_mat,0);
                drifts.a_F = max(model.d + obj.grids.a.matrix * (obj.p.r_a + obj.p.deathrate*obj.p.perfectannuities)...
                            + obj.p.directdeposit * obj.y_mat,0);
                drifts.b_B = min(model.s - model.d - aux.AdjustmentCost.cost(...
                	model.d, obj.grids.a.matrix, obj.p), 0);
                drifts.b_F = max(model.s - model.d - aux.AdjustmentCost.cost(...
                	model.d, obj.grids.a.matrix, obj.p), 0);
            elseif strcmp(obj.gridtype, 'HJB')
                drifts.a_B = min(model.d, 0)...
                	+ min(obj.grids.a.matrix * (obj.p.r_a + obj.p.deathrate * obj.p.perfectannuities)...
                	+ obj.p.directdeposit * obj.y_mat, 0);
                drifts.a_F = max(model.d, 0)...
                	+ max(obj.grids.a.matrix * (obj.p.r_a + obj.p.deathrate * obj.p.perfectannuities)...
                	+ obj.p.directdeposit * obj.y_mat, 0);
                drifts.b_B = min(-model.d - aux.AdjustmentCost.cost(model.d, obj.grids.a.matrix, obj.p), 0)...
                	+ min(model.s, 0);
                drifts.b_F = max(-model.d - aux.AdjustmentCost.cost(model.d, obj.grids.a.matrix, obj.p), 0)...
                	+ max(model.s, 0);    
            end
        end

    	%% --------------------------------------------------------------------
        % Liquid Asset Transitions
        % ---------------------------------------------------------------------
    	function A_liquid = compute_liquid_transitions(obj, bdriftB, bdriftF)
    		% Computes the sparse, square matrix for liquid asset transitions,
    		% adjusted for grid deltas.

    		lowdiag = - bdriftB ./ obj.grids.b.dB; 
            lowdiag(1,:,:,:) = 0;

            centerdiag = bdriftB ./ obj.grids.b.dB - bdriftF ./ obj.grids.b.dF; 
            centerdiag(1,:,:,:) = - bdriftF(1,:,:,:)./obj.grids.b.dF(1,:); 
            centerdiag(obj.nb,:,:,:) = bdriftB(obj.nb,:,:,:)./obj.grids.b.dB(obj.nb,:);

            updiag = bdriftF ./ obj.grids.b.dF; 
            updiag(obj.nb,:,:,:) = 0;

            A_liquid = obj.put_on_diags(lowdiag, centerdiag, updiag, 1);
        end

    	%% --------------------------------------------------------------------
        % Illiquid Asset Transitions
        % ---------------------------------------------------------------------
        function A_illiquid = compute_illiquid_transitions(obj, adriftB, adriftF)
        	% Computes the sparse, square matrix for illiquid asset transitions,
    		% adjusted for grid deltas.

    		lowdiag = -adriftB ./ obj.grids.a.dB; 
            lowdiag(:,1,:,:) = 0;

            centerdiag = adriftB ./ obj.grids.a.dB - adriftF ./ obj.grids.a.dF; 
            centerdiag(:,1,:,:) = -adriftF(:,1,:,:) ./ obj.grids.a.dF(:,1); 
            centerdiag(:,obj.na,:,:) = adriftB(:,obj.na,:,:) ./ obj.grids.a.dB(:,obj.na);

            updiag = adriftF ./ obj.grids.a.dF; 
            updiag(:,obj.na,:,:) = zeros(obj.nb,1,obj.nz,obj.ny);

            A_illiquid = obj.put_on_diags(lowdiag, centerdiag, updiag, 2);
        end

        %% --------------------------------------------------------------------
        % Computations for Returns Risk
        % ---------------------------------------------------------------------
        function perform_returns_risk_computations(obj)
            % Pre-computes terms needed for risky returns.
            %
            % Modified
            % --------
            % obj.asset_dB : deltaBackward for the grid of the risky asset
            %
            % obj.asset_dF : deltaForward for the grid of the risky asset
            %
            % obj.top : indicator for top of risky asset grid
            %
            % obj.risk_term : (asset * sigma_r) ^ 2
            %
            % obj.shift_dimension : dimension for offset of the risky asset on lower and upper diags
            %
            % obj.asset_dSum : deltaBackward + deltaForward

            obj.top = false(obj.nb, obj.na, obj.nz, obj.ny);
            if obj.p.OneAsset == 1
                obj.asset_dB = repmat(obj.grids.b.dB, [1 1 obj.nz obj.ny]);
                obj.asset_dF = repmat(obj.grids.b.dF, [1 1 obj.nz obj.ny]);

                obj.top(obj.nb, :, :, :) = true;

                obj.risk_term = (obj.grids.b.matrix * obj.p.sigma_r) .^ 2;

                obj.shift_dimension = 1;
            else
                obj.asset_dB = repmat(obj.grids.a.dB, [1 1 obj.nz obj.ny]);
                obj.asset_dF = repmat(obj.grids.a.dF, [1 1 obj.nz obj.ny]);

                obj.top(:, obj.na, :, :) = true;

                obj.risk_term = (obj.grids.a.matrix * obj.p.sigma_r) .^ 2;

                obj.shift_dimension = 2;
            end

            obj.asset_dF(obj.top) = obj.asset_dB(obj.top);
            obj.asset_dSum = obj.asset_dB + obj.asset_dF;
        end

        function Arisk_Vaa = compute_Vaa_terms(obj)
            % This function computes the (1/2) * (a * sigma_r) ^2 * Vaa component
            % of the A matrix.

            % Parameters
            % ----------
            % nb, na, nz, ny : grid lengths
            %
            % Returns
            % -------
            % Arisk_Vaa : the component of the A matrix s.t. Arisk_Vaa * V produces
            %   a second-order finite difference of V, after multiplying by the
            %   (1/2) * (a * sigma_r) ^ 2

            % (1/2) * (a * sigma_r) ^ 2 * V_{i+1} * (dF*(dB+dF)/2)
            V_i_plus_1_term = obj.risk_term ./ (obj.asset_dF .* obj.asset_dSum);

            % - (1/2) * (a * sigma_r) ^ 2 * V_i * (1/dB + 1/dF) / ((dB+dF)/2)
            V_i_term = - obj.risk_term .* (1./obj.asset_dB + 1./obj.asset_dF) ./ obj.asset_dSum;
            
            % (1/2) * (a * sigma_r) ^ 2 * V_{i-1} * (dB*(dB+dF)/2)
            V_i_minus_1_term = obj.risk_term ./ (obj.asset_dB .* obj.asset_dSum);

            % impose V' = 0 at the top of the grid
            V_i_term(obj.top) = - obj.risk_term(obj.top) ./ (obj.asset_dB(obj.top) .* obj.asset_dSum(obj.top));
            V_i_plus_1_term(obj.top) = 0;

            Arisk_Vaa = obj.put_on_diags(V_i_minus_1_term, V_i_term, V_i_plus_1_term, obj.shift_dimension);
        end

        function A = put_on_diags(obj, lower, middle, upper, shift_dim)
            % constructs a sparse matrix by putting column vectors on the diagonals

            % Parameters
            % ----------
            % lower : the vector to go on the lower diagonal
            %
            % middle : the vector to go on the middle diagonal
            %
            % upper : the vector to go on the upper diagonal
            %
            % shift : the offset of the lower and upper diagonals from the middle diagonal
            %
            % Returns
            % -------
            % A : a sparse matrix with lower, middle, and upper on the requested diagonals

            if shift_dim == 1
                upper = upper(:);
                upper = [0; upper(1:end-1)];
                lower = lower(:);
                lower = [lower(2:end); 0];
                A = spdiags(upper(:), 1, obj.n_states, obj.n_states) ...
                    + spdiags(middle(:), 0, obj.n_states, obj.n_states) ...
                    + spdiags(lower(:), -1, obj.n_states, obj.n_states);
            elseif shift_dim == 2
                upper = upper(:);
                upper = [zeros(obj.nb, 1); upper(1:end-obj.nb)];
                lower = lower(:);
                lower = [lower(obj.nb+1:end); zeros(obj.nb, 1)];
                A = spdiags(upper(:), obj.nb, obj.n_states, obj.n_states) ...
                    + spdiags(middle(:), 0, obj.n_states, obj.n_states) ...
                    + spdiags(lower(:), -obj.nb, obj.n_states, obj.n_states);                
            end
        end
    end
end

