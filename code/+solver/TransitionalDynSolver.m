classdef TransitionalDynSolver < handle
	% This superclass provides properties and methods 
	% for solving for policy functions given a future shock
	% by iterating backward on the dynamic HJB.
	%
	% MPCs out of news can be computed with this class
	% by iterating over the Feynman-Kac equation while
	% iterating over the dynamic HJB.
	%
	% This class allows for saving the policy functions
	% so they can be used to simulate MPCs out of news.

	properties (SetAccess = protected)
		p;
		income;
		grids;

        dim2;

		% intermediate values of V, A, policy fns, cumcon
		V;
		A;
		KFEint;
        
        shocks;
        savedTimesUntilShock;

		cumcon;
		ytrans_offdiag;

		% cumulative consumption for q1, q4 shock
		cum_con_q1 = cell(1,6);
		cum_con_q4 = cell(1,6);

		% results
		mpcs = struct();

		rho_diag;

	end

	methods
		function obj = TransitionalDynSolver(params,income,grids,shocks)

			obj.p = params;
			obj.income = income;
			obj.grids = grids;
            obj.shocks = shocks;

			obj.ytrans_offdiag = income.ytrans - diag(diag(income.ytrans));

			obj.mpcs = struct();
			for ishock = 1:6
				obj.mpcs(ishock).avg_1_quarterly = NaN;
				obj.mpcs(ishock).avg_4_quarterly = NaN(4,1);
				obj.mpcs(ishock).avg_4_annual = NaN;
            end
		end

		function solve(obj,KFE,pmf,cum_con_baseline)
            if obj.p.nz > 1
                warning('Transitional dynamics not coded for nz > 1')
                return
            end
            
            if obj.p.SimulateMPCS_news == 1
				savedTimesUntilShock = [4:-0.1:0.1 obj.p.delta_mpc];
				obj.savedTimesUntilShock = round(savedTimesUntilShock*40)/40;
				save([obj.p.tempdirec 'savedTimesUntilShock.mat'],'savedTimesUntilShock')
            end
            
            % loop over shocks
			for ishock = obj.shocks
                fprintf('    --- Shock = %f ---\n',obj.p.mpc_shocks(ishock))
				success = obj.getTerminalCondition(KFE,ishock);
				if ~success
					return
				end
				obj.iterateBackwards(ishock);

				if obj.p.ComputeMPCS_news == 1
					obj.computeMPCs(pmf,ishock,cum_con_baseline);
				end
            end
            fprintf('\n')
		end

		function success = getTerminalCondition(obj,KFE,ishock)
			% this method finds the value function the instant before
			% the shock is applied
			
            shock = obj.p.mpc_shocks(ishock);

            success = false;

			% Get the guess of terminal value function
			% V_{T+1}as V(b+shock,a,y)

			if strcmp(obj.dim2Identity,'a')
				dim2vec = obj.grids.a.vec;
				dim2mat = obj.grids.a.matrix;
			elseif strcmp(obj.dim2Identity,'c')
				dim2vec = obj.grids.c.vec;
				dim2mat = obj.grids.c.matrix;
			end

			if numel(p.rhos) > 1
		        rho_mat = reshape(obj.p.rhos,[1 1 numel(obj.p.rhos)]);
		    else
		        rho_mat = obj.p.rho;
		    end

			if (obj.p.ny > 1) && (obj.p.nz > 1)
				interp_grids = {obj.grids.b.vec,dim2vec,obj.grids.z.vec,obj.income.y.vec};
            elseif obj.p.ny > 1
                interp_grids = {obj.grids.b.vec,dim2vec,obj.income.y.vec};
            elseif obj.p.nz > 1
                interp_grids = {obj.grids.b.vec,dim2vec,obj.grids.z.vec};
			else
				interp_grids = {obj.grids.b.vec,dim2vec};
			end
			Vinterp = griddedInterpolant(interp_grids,squeeze(KFE.Vn),'linear');

			if (obj.p.ny > 1) && (obj.p.nz > 1)
				Vg_terminal = Vinterp(obj.grids.b.matrix(:)+shock,...
					dim2mat(:),obj.grids.z.matrix(:),obj.income.y.matrixKFE(:));
            elseif obj.p.ny > 1
                Vg_terminal = Vinterp(obj.grids.b.matrix(:)+shock,...
					dim2mat(:),obj.income.y.matrixKFE(:));
            elseif obj.p.nz > 1
                Vg_terminal = Vinterp(obj.grids.b.matrix(:)+shock,...
					dim2mat(:),obj.grids.z.matrix(:));
			else
				Vg_terminal = Vinterp(obj.grids.b.matrix(:)+shock,...
					dim2mat(:));
			end
			reshape_vec = [obj.p.nb_KFE,obj.dim2,obj.p.nz,obj.p.ny];
			Vg_terminal = reshape(Vg_terminal,reshape_vec);
			Vg_terminal_k = reshape(Vg_terminal,[],obj.p.ny);

			if strcmp(obj.dim2Identity,'c')
				KFE_terminal.c = obj.grids.c.matrix;
			end

			% iterate with implicit-explicit scheme to get V_terminal
			reshape_vec = [obj.p.nb_KFE,obj.dim2,obj.p.nz,obj.p.ny];
		    V_terminal = Vg_terminal;
		    V_terminal_k = reshape(V_terminal,[],obj.p.ny);
            
		    for ii = 1:5000
		    	if strcmp(obj.dim2Identity,'a')
			    	KFE_terminal = solver.find_policies(obj.p,obj.income,obj.grids,V_terminal);
                    A_terminal = solver.construct_trans_matrix(obj.p,obj.income,obj.grids,KFE_terminal,'KFE');
			    elseif strcmp(obj.dim2Identity,'c')
			    	KFE_terminal.h = solver.con_effort.find_policies(obj.p,obj.income,obj.grids,V_terminal);
	                KFE_terminal.s = (obj.p.r_b+obj.p.deathrate*obj.p.perfectannuities) * obj.grids.b.matrix...
	                    + (1-obj.p.wagetax)* obj.income.y.matrixKFE - obj.grids.c.matrix;

	                if p.SDU == 0
                    	KFE_terminal.u = aux.u_fn(obj.grids.c.matrix, obj.p.riskaver);
                    else
                    	KFE_terminal.u = aux.u_fn(obj.grids.c.matrix, obj.p.invies, rho_mat);
                    end
                    KFE_terminal.u = KFE_terminal.u - aux.con_effort.con_adj_cost(obj.p,KFE_terminal.h)...
	                    - aux.con_effort.penalty(obj.grids.b.matrix,obj.p.penalty1,obj.p.penalty2);

                    A_terminal = solver.con_effort.construct_trans_matrix(obj.p,obj.income,obj.grids,KFE_terminal);        
			    end
		    	u_k = reshape(KFE_terminal.u,[],obj.p.ny);

		    	V1_terminal_k = zeros(obj.p.nb_KFE*obj.dim2*obj.p.nz,obj.p.ny);
                deltaLarge = 1e-3;
		    	for k = 1:obj.p.ny
		    		ind1 = 1+obj.p.nb_KFE*obj.dim2*obj.p.nz*(k-1);
			    	ind2 = obj.p.nb_KFE*obj.dim2*obj.p.nz*k;
		    		Ak = A_terminal(ind1:ind2,ind1:ind2);
   
		    		indx_k = ~ismember(1:obj.p.ny,k);
		    		Vk_stacked 	= sum(repmat(obj.income.ytrans(k,indx_k),obj.p.nb_KFE*obj.dim2*obj.p.nz,1) ...
                            .* V_terminal_k(:,indx_k),2);

                    V1_terminal_k(:,k) = (deltaLarge* obj.rho_diag + (deltaLarge/obj.p.delta_mpc + 1 + deltaLarge*(obj.p.deathrate - obj.income.ytrans(k,k)))...
		        				*speye(obj.p.nb_KFE*obj.dim2*obj.p.nz) - deltaLarge* Ak)...
		                 	\ (deltaLarge* (u_k(:,k) + Vk_stacked)...
		                 		+ deltaLarge*Vg_terminal_k(:,k)/obj.p.delta_mpc + V_terminal_k(:,k));
		    	end

    			dst = max(abs(V1_terminal_k(:)-V_terminal(:)));
		        if (ii==1) || (mod(ii,10) == 0)
		            fprintf('\tFinding terminal value fn for mpc out of news, iter = %d, diff = %d\n',ii,dst)
		        end
		        
		        V_terminal_k = V1_terminal_k;
                V_terminal = reshape(V_terminal_k,reshape_vec);
		        
		        if dst < 1e-7
		        	fprintf('\tFound terminal value function after %i iterations\n',ii);
		        	success = true;
		            break
		        end
		    end

		    obj.V = V_terminal;
		    obj.KFEint = KFE_terminal;
            obj.A = A_terminal;
		end

		function iterateBackwards(obj,ishock)
			% iterate over the dynamic HJB
			
            shock = obj.p.mpc_shocks(ishock);
			obj.cumcon = zeros(obj.p.nb_KFE*obj.dim2*obj.income.ny*obj.p.nz,4);

			dim = obj.p.nb_KFE*obj.dim2*obj.income.ny*obj.p.nz;

			if strcmp(obj.dim2Identity,'c')
            	obj.KFEint.c = obj.grids.c.matrix;
            end

			for it = 4:-obj.p.delta_mpc:obj.p.delta_mpc
				timeUntilShock = obj.p.delta_mpc + 4 - it;
                timeUntilShock = round(timeUntilShock * 40) / 40;
				if mod(it,0.5) == 0
		            fprintf('\tUpdating cumulative consumption given news, quarter=%0.2f\n',it)
                end
                
                u_k = reshape(obj.KFEint.u,[],obj.p.ny);
                V_k = reshape(obj.V,[],obj.p.ny);
                V_k1 = zeros(size(V_k));
                for k = 1:obj.p.ny
                    ind1 = 1+obj.p.nb_KFE*obj.dim2*obj.p.nz*(k-1);
			    	ind2 = obj.p.nb_KFE*obj.dim2*obj.p.nz*k;
                    Ak = obj.A(ind1:ind2,ind1:ind2);
                    
                    indx_k = ~ismember(1:obj.income.ny,k);
		    		Vk_stacked 	= sum(repmat(obj.income.ytrans(k,indx_k),obj.p.nb_KFE*obj.dim2*obj.p.nz,1) ...
                            .* V_k(:,indx_k),2);
                        
                    V_k1(:,k) = (obj.rho_diag + (1/obj.p.delta_mpc + obj.p.deathrate - obj.income.ytrans(k,k))...
		        				* speye(obj.p.nb_KFE*obj.dim2*obj.p.nz) - Ak)...
		                 	\ (u_k(:,k) + Vk_stacked + V_k(:,k)/obj.p.delta_mpc);
                end
                obj.V = reshape(V_k1,[obj.p.nb_KFE,obj.dim2,obj.income.ny,obj.p.nz]);

		        % find policies a fraction of a period back
		        obj.update_policies();
		        
		        % find A matrix a fraction of a period back
		        obj.update_A_matrix();

			    if obj.p.ComputeMPCS_news == 1
				    % matrix divisor for Feynman-Kac
				    FKmat = cell(1,obj.income.ny);
			        for k = 1:obj.income.ny
			        	ind1 = 1+obj.p.nb_KFE*obj.dim2*obj.p.nz*(k-1);
			        	ind2 = obj.p.nb_KFE*obj.dim2*obj.p.nz*k;
			            FKmat{k} = speye(obj.p.nb_KFE*obj.dim2*obj.p.nz)*(...
			            					1/obj.p.delta_mpc + obj.p.deathrate - obj.income.ytrans(k,k))...
			                    	- obj.A(ind1:ind2,ind1:ind2);
			            FKmat{k} = inverse(FKmat{k});
			        end

			        for period = ceil(it):4
			        	obj.update_cum_con(period,FKmat);
	                end

			        if (it == 3 + obj.p.delta_mpc)
			        	obj.cum_con_q1{ishock} = obj.cumcon(:,4);
			        end
                end

			    if (obj.p.SimulateMPCS_news==1) && ismember(timeUntilShock,obj.savedTimesUntilShock)
			    	% save policy function
			    	index = find(obj.savedTimesUntilShock==timeUntilShock);
                    
                    obj.savePolicies(index,ishock);
			    	
                end
            end

            if obj.p.ComputeMPCS_news == 1
	            obj.cum_con_q4{ishock}(:,1) = obj.cumcon(:,1);
	            for period = 2:4
	                obj.cum_con_q4{ishock}(:,period) =...
	                	obj.cumcon(:,period) - obj.cumcon(:,period-1);
	            end
	        end
		end

		function update_cum_con(obj,period,FKmat)
            cumcon_t_k = reshape(obj.cumcon(:,period),[],obj.p.ny);
			for k = 1:obj.income.ny
                ytrans_cc_k = sum(obj.ytrans_offdiag(k,:) .* cumcon_t_k,2);

                deathin_cc_k = obj.get_death_inflows(cumcon_t_k,k);

                ind1 = 1+obj.p.nb_KFE*obj.dim2*obj.p.nz*(k-1);
                ind2 = obj.p.nb_KFE*obj.dim2*obj.p.nz*k;
                RHS = reshape(obj.KFEint.c(:,:,:,k),[],1) + ytrans_cc_k + deathin_cc_k ...
                        + cumcon_t_k(:,k)/obj.p.delta_mpc;
                obj.cumcon(ind1:ind2,period) = FKmat{k}*RHS;
            end
		end

		function computeMPCs(obj,pmf,ishock,cum_con_baseline)
            shock = obj.p.mpc_shocks(ishock);
            
            mpcs_1_quarterly = (obj.cum_con_q1{ishock} - cum_con_baseline(:,1)) / shock;
            obj.mpcs(ishock).avg_1_quarterly = mpcs_1_quarterly(:)' * pmf(:);

            mpcs_4_quarterly = (obj.cum_con_q4{ishock} - cum_con_baseline) / shock;
            obj.mpcs(ishock).avg_4_quarterly = mpcs_4_quarterly' * pmf(:);
            obj.mpcs(ishock).avg_4_annual = sum(obj.mpcs(ishock).avg_4_quarterly);
		end
	end

end