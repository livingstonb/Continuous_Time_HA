classdef TransitionMatrixConstructorSDU < HACTLib.computation.TransitionMatrixConstructor

	methods
		function [A, stationary] = construct(obj, model, V)
			[A, stationary, drifts] = construct@HACTLib.computation.TransitionMatrixConstructor(obj, model, V);

			if obj.returns_risk
            	if obj.p.OneAsset
            		A = A + obj.compute_Va_terms(drifts.b_B, drifts.b_F, V);
	            else
	            	A = A + obj.compute_Va_terms(drifts.a_B, drifts.a_F, V);
	            end
            end 
		end
	end

	methods (Access=protected)
		function [Arisk_Va, stationary] = compute_Va_terms(obj, driftB, driftF, V)
            % Computes the (1/2) * (a * sigma_r) ^ 2 * Va ^2 ... term when
            % returns are risky.
            %
            % Parameters
            % ----------
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

            assert(obj.p.SDU, "This term only applies to SDU")

            [V1B, V1F] = obj.first_diffs(V);
            V1 = (driftB < 0) .* V1B + (driftF > 0) .* V1F;

            % Now add zeta(a) * Va^2 term.
            adj_term = make_SDU_adjustment(obj, V, V1);

            updiag = zeros(obj.shape);
            lowdiag = zeros(obj.shape);

            lowdiag(driftB < 0) = - adj_term(driftB < 0) ./ obj.asset_dB(driftB < 0);
            updiag(driftF > 0) = adj_term(driftF > 0) ./ obj.asset_dF(driftF > 0);

            if obj.p.OneAsset
                lowdiag(1, :, :, :) = 0;
                updiag(obj.nb, :, :, :) = 0;
            else
                lowdiag(:, 1, :, :) = 0;
                updiag(:, obj.na, :, :) = 0;
            end

            centdiag = - lowdiag - updiag;
            Arisk_Va = HACTLib.aux.sparse_diags(...
                [lowdiag(:), centdiag(:), updiag(:)],...
                obj.offsets_for_rr);

            % Some states may have no drift, we will need to deal with them
            % separately by adding directly to the right_hand_side of the HJB.
            stationary = (driftB >= 0) & (driftF <= 0);
        end

        function [V1B, V1F] = first_diffs(obj, V)
        	V1B = zeros(obj.shape);
            V1F = zeros(obj.shape);
            
            if obj.p.OneAsset
                V1B(2:obj.nb,:,:,:) = (V(2:obj.nb,:,:,:) - V(1:obj.nb-1,:,:,:))...
                	./ obj.grids.b.dB(2:obj.nb,:);
                V1F(1:obj.nb-1,:,:,:) = (V(2:obj.nb,:,:,:) - V(1:obj.nb-1,:,:,:))...
                	./ obj.grids.b.dF(1:obj.nb-1,:);
            else
                V1B(:,2:obj.na,:,:) = (V(:,2:obj.na,:,:) - V(:,1:obj.na-1,:,:))...
                	./ obj.grids.a.dB(:,2:obj.na);
                V1F(:,1:obj.na-1,:,:) = (V(:,2:obj.na,:,:) - V(:,1:obj.na-1,:,:))...
                	./ obj.grids.a.dF(:,1:obj.na-1);
            end
        end

        function SDU_adj = make_SDU_adjustment(obj, V, V1)
        	if obj.p.invies == 1
                SDU_adj = obj.risk_term .* V1 * (1 - obj.p.riskaver);
            else
                SDU_adj = obj.risk_term .* V1 ./...
                	V * (obj.p.invies - obj.p.riskaver) / (1 - obj.p.invies);
            end
        end
	end
end