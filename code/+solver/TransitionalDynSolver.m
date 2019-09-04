classdef TransitionalDynSolver < handle

	properties (SetAccess = protected)
		p;
		income;
		grids;

        dim2;

		% intermediate values of V, A, policy fns, cumcon
		V;
		A;
		KFEint;
        
        shocks = [5];
        savedTimesUntilShock;

		cumcon;
		ytrans_offdiag;

		% cumulative consumption for q1, q4 shock
		cum_con_q1 = cell(1,6);
		cum_con_q4 = cell(1,6);

		% results
		mpcs = struct();

	end

	methods
		function obj = TransitionalDynSolver(params,income,grids)
			obj.p = params;
			obj.income = income;
			obj.grids = grids;

			obj.ytrans_offdiag = income.ytrans - diag(diag(income.ytrans));

			obj.mpcs = struct();
			for ishock = 1:6
				obj.mpcs(ishock).avg_1_t = NaN;
				obj.mpcs(ishock).avg_4_t = NaN(4,1);
			end
		end

		function solve(obj,KFE,pmf,cum_con_baseline)
            if obj.p.SimulateMPCS_news == 1
				savedTimesUntilShock = [4:-0.2:0.2 obj.p.delta_mpc];
				obj.savedTimesUntilShock = round(savedTimesUntilShock*40)/40;
				save([obj.p.tempdirec 'savedTimesUntilShock.mat'],'savedTimesUntilShock')
            end
            
            % loop over shocks
			for ishock = obj.shocks
                fprintf('    --- Shock = %f ---\n',obj.p.mpc_shocks(ishock))
				obj.getTerminalCondition(KFE,ishock);
				obj.iterateBackwards(ishock);

				if obj.p.ComputeMPCS == 1
					obj.computeMPCs(pmf,ishock,cum_con_baseline);
				end
            end
            fprintf('\n')
		end

		function getTerminalCondition(obj,KFE,ishock)
			% this method finds the value function the instant before
			% the shock is applied
			
            shock = obj.p.mpc_shocks(ishock);

			% Get the guess of terminal value function
			% V_{T+1}as V(b+shock,a,y)

			if strcmp(obj.dim2Identity,'a')
				dim2vec = obj.grids.a.vec;
				dim2mat = obj.grids.a.matrix;
			elseif strcmp(obj.dim2Identity,'c')
				dim2vec = obj.grids.c.vec;
				dim2mat = obj.grids.c.matrix;
			end

			if obj.p.ny > 1
				interp_grids = {obj.grids.b.vec,dim2vec,obj.income.y.vec};
			else
				interp_grids = {obj.grids.b.vec,dim2vec};
			end
			Vinterp = griddedInterpolant(interp_grids,KFE.Vn,'linear');

			if obj.p.ny > 1
				Vg_terminal = Vinterp(obj.grids.b.matrix(:)+shock,...
					dim2mat(:),obj.income.y.matrixKFE(:));
			else
				Vg_terminal = Vinterp(obj.grids.b.matrix(:)+shock,...
					dim2mat(:));
			end
			reshape_vec = [obj.p.nb_KFE,obj.dim2,obj.p.ny,obj.p.nz];
			Vg_terminal = reshape(Vg_terminal,reshape_vec);
			Vg_terminal = permute(Vg_terminal,[1 2 4 3]);
			Vg_terminal_k = reshape(Vg_terminal,[],obj.p.ny);

			% iterate with implicit-explicit scheme to get V_terminal
			reshape_vec_switched = [obj.p.nb_KFE,obj.dim2,obj.p.nz,obj.p.ny];
		    V_terminal = Vg_terminal;
		    V_terminal_k = reshape(Vg_terminal,[],obj.p.ny);
		    for ii = 1:5000
		    	V_normaldim = permute(V_terminal,[1 2 4 3]);
		    	KFE_terminal = solver.two_asset.find_policies(obj.p,obj.income,obj.grids,V_normaldim);
		    	A_terminal = solver.two_asset.construct_trans_matrix(obj.p,obj.income,obj.grids,KFE_terminal,'KFE');
		    	u_switched = permute(KFE_terminal.u,[1 2 4 3]);
		    	u_switched_k = reshape(u_switched,[],obj.p.ny);

		    	V1_terminal_k = zeros(obj.p.nb_KFE*obj.dim2*obj.p.nz,obj.p.ny);
		    	for k = 1:obj.p.ny
		    		ind1 = 1+obj.p.nb_KFE*obj.dim2*obj.p.nz*(k-1);
			    	ind2 = obj.p.nb_KFE*obj.dim2*obj.p.nz*k;
		    		Ak = A_terminal(ind1:ind2,ind1:ind2);

		    		indx_k = ~ismember(1:obj.p.ny,k);
		    		Vk_stacked 	= sum(repmat(obj.income.ytrans(k,indx_k),obj.p.nb_KFE*obj.dim2*obj.p.nz,1) ...
                            .* V_terminal_k(:,indx_k),2);

		    		V1_terminal_k(:,k) = ((1/obj.p.delta_mpc + 1/(1e-4) + obj.p.rho + obj.p.deathrate - obj.income.ytrans(k,k))...
		        				*speye(obj.p.nb_KFE*obj.dim2*obj.p.nz) - Ak)...
		                 	\ (u_switched_k(:,k) + Vk_stacked...
		                 		+ Vg_terminal_k(:,k)/(1e-4) + V_terminal_k(:,k)/obj.p.delta_mpc);
		    	end

    			dst = max(abs(V1_terminal_k(:)-V_terminal(:)));
		        if (ii==1) || (mod(ii,10) == 0)
		            fprintf('\tFinding terminal value fn for mpc out of news, iter = %d, diff = %d\n',ii,dst)
		        end
		        
		        V_terminal_k = V1_terminal_k;
		    	V_terminal = reshape(V1_terminal_k,reshape_vec_switched);
		        
		        if dst < 1e-10
		        	fprintf('\tFound terminal value function after %i iterations\n',ii)
		            break
		        end
		    end

		    obj.V = V_terminal;
		    obj.KFEint = KFE_terminal;
            obj.A = A_terminal;
		end

		function iterateBackwards(obj,ishock)
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
                
                u_switch = permute(obj.KFEint.u,[1 2 4 3]);
                u_switch_k = reshape(u_switch,[],obj.p.ny);
                V_k = reshape(obj.V,[],obj.p.ny);
                V_k1 = zeros(size(V_k));
                for k = 1:obj.p.ny
                    ind1 = 1+obj.p.nb_KFE*obj.dim2*(k-1);
			    	ind2 = obj.p.nb_KFE*obj.dim2*k;
                    Ak = obj.A(ind1:ind2,ind1:ind2);
                    
                    indx_k = ~ismember(1:obj.income.ny,k);
		    		Vk_stacked 	= sum(repmat(obj.income.ytrans(k,indx_k),obj.p.nb_KFE*obj.dim2,1) ...
                            .* V_k(:,indx_k),2);
                        
                    V_k1(:,k) = ((1/obj.p.delta_mpc + obj.p.rho + obj.p.deathrate - obj.income.ytrans(k,k))...
		        				* speye(obj.p.nb_KFE*obj.dim2) - Ak)...
		                 	\ (u_switch_k(:,k) + Vk_stacked...
		                 		+ V_k(:,k)/obj.p.delta_mpc);
                end

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
			for k = 1:obj.income.ny
                cumcon_t = obj.cumcon(:,period);
                ytrans_cc_k = sum(obj.ytrans_offdiag(k,:) ...
                    .* reshape(cumcon_t,[],obj.income.ny),2);

                deathin_cc_k = obj.get_death_inflows(cumcon_t,k);

                ind1 = 1+obj.p.nb_KFE*obj.dim2*(k-1);
                ind2 = obj.p.nb_KFE*obj.dim2*k;
                RHS = reshape(obj.KFEint.c(:,:,k,:),[],1) + ytrans_cc_k + deathin_cc_k ...
                        + cumcon_t(ind1:ind2)/obj.p.delta_mpc;
                obj.cumcon(ind1:ind2,period) = FKmat{k}*RHS;
            end
		end

		function computeMPCs(obj,pmf,ishock,cum_con_baseline)
            shock = obj.p.mpc_shocks(ishock);
            
            mpcs_1_t = (obj.cum_con_q1{ishock} - cum_con_baseline(:,1)) / shock;
            obj.mpcs(ishock).avg_1_t = mpcs_1_t(:)' * pmf(:);

            mpcs_4_t = (obj.cum_con_q4{ishock} - cum_con_baseline) / shock;
            obj.mpcs(ishock).avg_4_t = mpcs_4_t' * pmf(:);
		end
	end

end