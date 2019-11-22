classdef A_Matrix_Constructor < handle
    % This class constructs the A matrix for the HJB or KFE
    % First, instantiate the class with the required arguments,
    % then call the construct method with the policy functions
    % and value function as arguments.
    properties (SetAccess = protected)
        % grid sizes
        nb;
        na;
        nz;
        ny;
        dim;
       
        % grids and parameters
        grids;
        p;

        % returns risk is on or off
        returns_risk;

        % grid deltas
        asset_dB;
        asset_dF;
        asset_dSum; % (deltaB + deltaF)

        % boolean mask for top of the b (1-asset) or a (2-asset) grid
        top;

        % pre-computed term when returns are risky
        risk_term;

        % offset of terms on lower/upper diagonals for risky asset
        % 1 for liquid, 2 for illiquid
        shift_dimension;

        % income values
        y_mat;

        % 'HJB' or 'KFE' grid
        modeltype;
    end

    methods
        %% --------------------------------------------------------------------
        % CLASS CONSTRUCTOR
        % ---------------------------------------------------------------------
        function obj = A_Matrix_Constructor(p, income, grids, modeltype, returns_risk)
            % instantiate the constructor

            % Parameters
            % ----------
            % p : a Params object which contains the parameters of the model
            %
            % income : an Income object
            %
            % grids : a Grid object which contains the grids
            %
            % modeltype : a string ('HJB' or 'KFE') indicating which grids are
            %   to be used
            %
            % returns_risk : a flag that determines whether returns risk is to
            %   be included in the A matrix
            %
            % Returns
            % -------
            % an A_Matrix_Constructor object

            obj.nb = numel(grids.b.vec);
            obj.na = numel(grids.a.vec);
            obj.nz = p.nz;
            obj.ny = numel(income.y.vec);
            obj.dim = obj.nb * obj.na * obj.nz * obj.ny;

            obj.grids = grids;
            
            obj.p = p;

            obj.modeltype = modeltype;

            obj.returns_risk = returns_risk;

            if strcmp(modeltype,'KFE')
                obj.y_mat = income.y.matrixKFE;
            elseif strcmp(modeltype,'HJB')
                obj.y_mat = income.y.matrix;
            end

            if p.sigma_r > 0
                obj.perform_returns_risk_computations();
            end
        end
        
        %% --------------------------------------------------------------------
        % COMPUTATIONS FOR RETURNS RISK
        % ---------------------------------------------------------------------
        function perform_returns_risk_computations(obj)
            % pre-computes terms needed for risky returns

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

        %% --------------------------------------------------------------------
        % CONSTRUCT THE A MATRIX
        % ---------------------------------------------------------------------
        function [A, stationary] = construct(obj, model, V)
            % constructs the transition matrix

            % Parameters
            % ----------
            % model : a structure containing the policy functions. saving (s)
            %   and deposits (d) are used here
            % V : the value function, shape (nb, na, nz, ny)
            %
            % Returns
            % -------
            % A : the sparse transition matrix of shape (nb*na*nz*ny, nb*na*nz*ny)
            %
            % stationary : a boolean mask indicating states where drift in the
            %   risky asset was neither backward nor forward. risk computations
            %   for these states will be included outside the A matrix

            % grids
            nb = obj.nb;
            na = obj.na;
            nz = obj.nz;
            ny = obj.ny;
            dim = obj.dim;

            % policy functions
            s = model.s;
            d = model.d;

            % states where drift in risky asset is neither forward nor backward
            stationary = [];

            %% ----------------------------------------------
            % INITIALIZE
            % -----------------------------------------------
            chi = NaN(nb,na,nz,ny);
            yy = NaN(nb,na,nz,ny);
            zetta = NaN(nb,na,nz,ny);

            X = NaN(nb,na,nz,ny);
            Y = NaN(nb,na,nz,ny);
            Z = NaN(nb,na,nz,ny);

            %% --------------------------------------------------------------------
            % COMPUTE ASSET DRIFTS
            % ---------------------------------------------------------------------
            if strcmp(obj.modeltype,'KFE')
                adriftB = min(d + obj.grids.a.matrix * (obj.p.r_a + obj.p.deathrate*obj.p.perfectannuities)...
                                                         + obj.p.directdeposit * obj.y_mat,0);
                adriftF = max(d + obj.grids.a.matrix * (obj.p.r_a + obj.p.deathrate*obj.p.perfectannuities)...
                                                         + obj.p.directdeposit * obj.y_mat,0);

                bdriftB = min(s - d - aux.AdjustmentCost.cost(d,obj.grids.a.matrix,obj.p),0);
                bdriftF = max(s - d - aux.AdjustmentCost.cost(d,obj.grids.a.matrix,obj.p),0);
            elseif strcmp(obj.modeltype,'HJB')
                adriftB = min(d,0) + min(obj.grids.a.matrix * (obj.p.r_a + obj.p.deathrate*obj.p.perfectannuities) + obj.p.directdeposit * obj.y_mat,0);
                adriftF = max(d,0) + max(obj.grids.a.matrix * (obj.p.r_a + obj.p.deathrate*obj.p.perfectannuities) + obj.p.directdeposit * obj.y_mat,0);

                bdriftB = min(-d - aux.AdjustmentCost.cost(d,obj.grids.a.matrix,obj.p),0) + min(s,0);
                bdriftF = max(-d - aux.AdjustmentCost.cost(d,obj.grids.a.matrix,obj.p),0) + max(s,0);    
            end

            %% --------------------------------------------------------------------
            % ILLIQUID ASSET TRANSITIONS
            % ---------------------------------------------------------------------
            chi(:,2:na,:,:) = -adriftB(:,2:na,:,:)./obj.grids.a.dB(:,2:na); 
            chi(:,1,:,:) = zeros(nb,1,nz,ny);
            yy(:,2:na-1,:,:) = adriftB(:,2:na-1,:,:)./obj.grids.a.dB(:,2:na-1) ...
                - adriftF(:,2:na-1,:,:)./obj.grids.a.dF(:,2:na-1); 
            yy(:,1,:,:) = - adriftF(:,1,:,:)./obj.grids.a.dF(:,1); 
            yy(:,na,:,:) = adriftB(:,na,:,:)./obj.grids.a.dB(:,na);
            zetta(:,1:na-1,:,:) = adriftF(:,1:na-1,:,:)./obj.grids.a.dF(:,1:na-1); 
            zetta(:,na,:,:) = zeros(nb,1,nz,ny);

            A = obj.put_on_diags(chi, yy, zetta, 2);

            %% --------------------------------------------------------------------
            % LIQUID ASSET TRANSITIONS
            % ---------------------------------------------------------------------
            X(2:nb,:,:,:) = - bdriftB(2:nb,:,:,:)./obj.grids.b.dB(2:nb,:); 
            X(1,:,:,:) = zeros(1,na,nz,ny);
            Y(2:nb-1,:,:,:) = bdriftB(2:nb-1,:,:,:)./obj.grids.b.dB(2:nb-1,:) ...
                - bdriftF(2:nb-1,:,:,:)./obj.grids.b.dF(2:nb-1,:); 
            Y(1,:,:,:) = - bdriftF(1,:,:,:)./obj.grids.b.dF(1,:); 
            Y(nb,:,:,:) = bdriftB(nb,:,:,:)./obj.grids.b.dB(nb,:);
            Z(1:nb-1,:,:,:) = bdriftF(1:nb-1,:,:,:)./obj.grids.b.dF(1:nb-1,:); 
            Z(nb,:,:,:) = zeros(1,na,nz,ny);

            A = A + obj.put_on_diags(X, Y, Z, 1);
            
            %% --------------------------------------------------------------------
            % RATE OF RETURN RISK
            % ---------------------------------------------------------------------
            if (obj.p.OneAsset == 1) && (obj.returns_risk == 1)

                % Second derivative component, Vaa
                A = A + obj.compute_V_aa_terms(nb, na, nz, ny);

                if obj.p.SDU == 1
                    % First derivative component, Va
                    A = A + obj.compute_V_a_terms(nb, na, nz, ny, V, adriftB, adriftF);
                end

            elseif (obj.p.OneAsset == 0) && (obj.returns_risk == 1)

                % Second derivative component, Vaa
                A = A + obj.compute_V_aa_terms(nb, na, nz, ny);

                if obj.p.SDU == 1
                    % First derivative component, Va
                    [Arisk_Va, stationary] = obj.compute_V_a_terms(nb, na, nz, ny, V, adriftB, adriftF);
                    A = A + Arisk_Va;
                end     
            end
        end

        function Arisk_Vaa = compute_V_aa_terms(obj, nb, na, nz, ny)
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

        function [Arisk_Va, stationary] = compute_V_a_terms(obj, nb, na,...
            nz, ny, V, driftB, driftF)
            % computes the (1/2) * (a * sigma_r) ^ 2 * Va ^2 ... term when returns are risky

            % Parameters
            % ----------
            % nb, na, nz, ny : grid lengths
            %
            % V : the value function, of shape (nb, nz, nz, ny)
            %
            % driftB : the backward drift of the risky asset, where negative
            %
            % driftF : the forward drift of the risky asset, where positive
            %
            % Returns
            % -------
            % Arisk_Va : after multiplying by V, this is the component of the A matrix that
            %   produces the Va ^ 2 term
            %
            % stationary : a boolean mask indicating states where drift in the
            %   risky asset was neither backward nor forward. risk computations
            %   for these states will be included outside the A matrix

            assert(obj.p.SDU == 1, "This term only applies to SDU.")

            V1B = zeros(nb, na, nz, ny);
            V1F = zeros(nb, na, nz, ny);
            
            if obj.p.OneAsset == 1
                V1B(2:nb,:,:,:) = (V(2:nb,:,:,:) - V(1:nb-1,:,:,:)) ./ obj.grids.b.dB(2:nb, :);
                V1F(1:nb-1,:,:,:) = (V(2:nb,:,:,:) - V(1:nb-1,:,:,:)) ./ obj.grids.b.dF(1:nb-1, :);
            else
                V1B(:,2:na,:,:) = (V(:,2:na,:,:) - V(:,1:na-1,:,:)) ./ obj.grids.a.dB(:,2:na);
                V1F(:,1:na-1,:,:) = (V(:,2:na,:,:) - V(:,1:na-1,:,:)) ./ obj.grids.a.dF(:,1:na-1);
            end

            V1 = (driftB < 0) .* V1B + (driftF > 0) .* V1F;

            % some states may have no drift, will need to deal with them separately
            % by adding to the RHS of HJB
            stationary = (driftB >= 0) & (driftF <= 0);

            % now add zeta(a) * Va^2 term
            if obj.p.invies == 1
                adj_term = obj.risk_term .* V1 * (1 - obj.p.riskaver);
            else
                adj_term = obj.risk_term .* V1 ./ V * (obj.p.invies - obj.p.riskaver) / (1 - obj.p.invies);
            end

            updiag = zeros(nb, na, nz, ny);
            lowdiag = zeros(nb, na, nz, ny);

            lowdiag(driftB < 0) = - adj_term(driftB < 0) ./ obj.asset_dB(driftB < 0);
            updiag(driftF > 0) = adj_term(driftF > 0) ./ obj.asset_dF(driftF > 0);

            if obj.p.OneAsset == 1
                lowdiag(1, :, :, :) = 0;
                updiag(nb, :, :, :) = 0;
            else
                lowdiag(:, 1, :, :) = 0;
                updiag(:, na, :, :) = 0;
            end
            centdiag = - lowdiag - updiag;

            Arisk_Va = obj.put_on_diags(lowdiag, centdiag, updiag, obj.shift_dimension);
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
                A = spdiags(upper(:), 1, obj.dim, obj.dim) ...
                    + spdiags(middle(:), 0, obj.dim, obj.dim) ...
                    + spdiags(lower(:), -1, obj.dim, obj.dim);
            elseif shift_dim == 2
                upper = upper(:);
                upper = [zeros(obj.nb, 1); upper(1:end-obj.nb)];
                lower = lower(:);
                lower = [lower(obj.nb+1:end); zeros(obj.nb, 1)];
                A = spdiags(upper(:), obj.nb, obj.dim, obj.dim) ...
                    + spdiags(middle(:), 0, obj.dim, obj.dim) ...
                    + spdiags(lower(:), -obj.nb, obj.dim, obj.dim);                
            end
        end
    end
end

