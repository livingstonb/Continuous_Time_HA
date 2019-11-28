classdef TransitionalDynSolver < handle
	% this class provides properties and methods 
	% for solving for policy functions given a future shock
	% by iterating backward on the dynamic HJB.
	%
	% MPCs out of news can be computed with this class
	% by iterating over the Feynman-Kac equation while
	% iterating over the dynamic HJB.
	%
	% this class allows for saving the policy functions
	% so they can be used to simulate MPCs out of news as well

	properties (SetAccess = protected)

		% parameters, grids of the model
		p;
		income;
		grids;

		% intermediate values of V, A, policy fns, cumcon
		V;
		A_constructor_HJB;
		A_constructor_FK;
		A_HJB;
		A_FK;
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
		function obj = TransitionalDynSolver(params, income, grids, shocks)
			% class constructor

			% Parameters
			% ----------
			% params : a Params object
			%
			% income : an Income object
			%
			% grids : a Grid object
			%
			% shocks : a vector of shock indices, in reference to the mpc shocks
			%	in 'params'
			%
			% Returns
			% -------
			% obj : a TransitionalDynSolver object

			% store references to important model objects
			obj.p = params;
			obj.income = income;
			obj.grids = grids;
            obj.shocks = shocks;

			obj.ytrans_offdiag = income.ytrans - diag(diag(income.ytrans));

			% initialize the results
			obj.mpcs = struct();
			for ishock = 1:6
				obj.mpcs(ishock).avg_1_quarterly = NaN;
				obj.mpcs(ishock).avg_4_quarterly = NaN(4,1);
				obj.mpcs(ishock).avg_4_annual = NaN;
            end

            % initialize constructors of the A matrix
            returns_risk = (obj.p.sigma_r > 0);
		    obj.A_constructor_HJB = HACTLib.computation.TransitionMatrixConstructor(...
                obj.p, obj.income, obj.grids, 'KFE', returns_risk);

		    % if returns are risky and returns risk is used in the KFE,
		    % this A constructor will be used instead:
		    if returns_risk && (obj.p.retrisk_KFE == 0)
		    	obj.A_constructor_FK = solver.A_Matrix_Constructor(obj.p, obj.income, obj.grids, 'KFE', false);
		    end
		end

		function solve(obj, KFE, pmf, cum_con_baseline)
			% calls various class methods to solve the model backward
			% and compute MPCs

			% Parameters
			% ----------
			% KFE : a structure containing baseline policy functions and
			%	the value function on the KFE grid
			%
			% pmf : the equilibrium distribution over the KFE grid, of
			%	shape (nb_KFE, na_KFE, nz, ny)
			%
			% cum_con_baseline : cumulative consumption for each period
			%	in the baseline, used to compute MPCs
			%
			% Directly Modifies
			% -----------------
			% obj.savedTimesUntilShock : a vector indicating the time until
			%	the shock at a given iteration, used to find the appropriate
			%	policy function given the current subperiod; only modified/used
			%	if SimulateMPCS_news is turned on
			%
			% Note
			% ----
			% This method calls various other class methods and so indirectly
			% modifies other class properties.

            if obj.p.nz > 1
                warning('Transitional dynamics not coded for nz > 1')
                return
            elseif (obj.p.SDU == 1) || (obj.p.sigma_r > 0)
            	warning('Transitional dynamics not coded for SDU or returns risk')
                return
            end
            
            if obj.p.SimulateMPCS_news == 1
				savedTimesUntilShock = [4:-0.1:0.1 obj.p.delta_mpc];
				obj.savedTimesUntilShock = round(savedTimesUntilShock*40)/40;
				save([obj.p.tempdirec 'savedTimesUntilShock.mat'],'savedTimesUntilShock')
            end

            if numel(obj.p.rhos) > 1
		        rhocol = kron(obj.p.rhos',ones(obj.p.nb_KFE*obj.p.na_KFE,1));
		        obj.rho_diag = spdiags(rhocol,0,obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz, obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz);
		    else
		        obj.rho_diag = obj.p.rho * speye(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz);
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

		function success = getTerminalCondition(obj, KFE, ishock)
			% this method finds the value function the instant before
			% the shock is applied

			% Parameters
			% ----------
			% KFE : structure containing the policy functions and the value
			%	function on the KFE grid
			%
			% ishock : the index of the shock, in reference to the vector
			%	of shock sizes in the Params object used to instantiate this
			%	class
			%
			% Returns
			% -------
			% success : a flag taking the value of true if the terminal
			%	condition was found, false otherwise
			%
			% Modifies
			% --------
			% obj.V : the value function given a future shock; this method
			%	initializes the value function with the terminal value
			%	function, of shape (nb_KFE, na_KFE, nz, ny)
			%
			% obj.KFEint : a structure containing the policy functions on
			%	the KFE grids; this method initializes these policy functions
			%	with the terminal policy functions
			%
			% obj.A : the transition matrix over the KFE grids; this method
			%	initializes A with the terminal transition matrix

            shock = obj.p.mpc_shocks(ishock);

            success = false;

			% Get the guess of terminal value function
			% V_{T+1}as V(b+shock,a,y)

			if numel(obj.p.rhos) > 1
		        rho_mat = reshape(obj.p.rhos,[1 1 numel(obj.p.rhos)]);
		    else
		        rho_mat = obj.p.rho;
		    end

			Vg_terminal = obj.get_value_at_instant_of_shock(KFE, shock);

			reshape_vec = [obj.p.nb_KFE,obj.p.na_KFE,obj.p.nz,obj.p.ny];
			Vg_terminal = reshape(Vg_terminal,reshape_vec);
			Vg_terminal_k = reshape(Vg_terminal,[],obj.p.ny);

			% iterate with implicit-explicit scheme to get V_terminal
			reshape_vec = [obj.p.nb_KFE,obj.p.na_KFE,obj.p.nz,obj.p.ny];
		    V_terminal = Vg_terminal;
		    V_terminal_k = reshape(V_terminal,[],obj.p.ny);
            
		    for ii = 1:5000
		    	KFE_terminal = solver.find_policies(obj.p,obj.income,obj.grids,V_terminal);
		    	A_terminal = obj.A_constructor_HJB.construct(KFE_terminal, V_terminal);

		    	u_k = reshape(KFE_terminal.u,[],obj.p.ny);

		    	V1_terminal_k = zeros(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz,obj.p.ny);
                deltaLarge = 1e-3;
		    	for k = 1:obj.p.ny
		    		ind1 = 1+obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*(k-1);
			    	ind2 = obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*k;
		    		Ak = A_terminal(ind1:ind2, ind1:ind2);
                    
                    ez_adj = obj.income.SDU_income_risk_adjustment(obj.p, V_terminal);
   
		    		indx_k = ~ismember(1:obj.p.ny,k);
		    		if obj.p.SDU == 0
			    		Vk_stacked 	= sum(repmat(obj.income.ytrans(k,indx_k),obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz,1) ...
	                            .* V_terminal_k(:,indx_k),2);
			    	else
			    		Vk_stacked = sum(squeeze(ez_adj(:,k,indx_k)) .* V_terminal_k(:,indx_k),2);
			    	end

			    	if obj.p.SDU == 0
	                    V1_terminal_k(:,k) = (deltaLarge* obj.rho_diag + (deltaLarge/obj.p.delta_mpc + 1 + deltaLarge*(obj.p.deathrate - obj.income.ytrans(k,k)))...
			        				*speye(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz) - deltaLarge* Ak)...
			                 	\ (deltaLarge* (u_k(:,k) + Vk_stacked)...
			                 		+ deltaLarge*Vg_terminal_k(:,k)/obj.p.delta_mpc + V_terminal_k(:,k));
			        else
			        	V1_terminal_k(:,k) = ...
			        		(deltaLarge* obj.rho_diag + (deltaLarge/obj.p.delta_mpc + 1 + ...
			        			deltaLarge*obj.p.deathrate.*speye(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz)) - deltaLarge* Ak...
			        			- deltaLarge * spdiags(ez_adj(:,k,k), 0, obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz, obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz))...
			                 	\ (deltaLarge* (u_k(:,k) + Vk_stacked)...
			                 		+ deltaLarge*Vg_terminal_k(:,k)/obj.p.delta_mpc + V_terminal_k(:,k));
		        	end
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
            obj.A_HJB = A_terminal;

            if (obj.p.sigma_r > 0) && (obj.p.retrisk_KFE == 0)
				obj.A_FK = obj.A_constructor_FK.construct(obj.KFEint);
			end
		end

		function V_at_shock = get_value_at_instant_of_shock(obj, KFE, shock)
			% V_at_shock is needed to find the value function an instant
			% before the shock

			if (obj.p.ny > 1) && (obj.p.nz > 1)
				interp_grids = {obj.grids.b.vec,obj.grids.a.vec,obj.grids.z.vec,obj.income.y.vec};
            elseif obj.p.ny > 1
                interp_grids = {obj.grids.b.vec,obj.grids.a.vec,obj.income.y.vec};
            elseif obj.p.nz > 1
                interp_grids = {obj.grids.b.vec,obj.grids.a.vec,obj.grids.z.vec};
			else
				interp_grids = {obj.grids.b.vec,obj.grids.a.vec};
			end
			Vinterp = griddedInterpolant(interp_grids,squeeze(KFE.Vn),'linear');

			if (obj.p.ny > 1) && (obj.p.nz > 1)
				V_at_shock = Vinterp(obj.grids.b.matrix(:)+shock,...
					obj.grids.a.matrix(:),obj.grids.z.matrix(:),obj.income.y.matrixKFE(:));
            elseif obj.p.ny > 1
                V_at_shock = Vinterp(obj.grids.b.matrix(:)+shock,...
					obj.grids.a.matrix(:),obj.income.y.matrixKFE(:));
            elseif obj.p.nz > 1
                V_at_shock = Vinterp(obj.grids.b.matrix(:)+shock,...
					obj.grids.a.matrix(:),obj.grids.z.matrix(:));
			else
				V_at_shock = Vinterp(obj.grids.b.matrix(:)+shock,...
					obj.grids.a.matrix(:));
			end
		end

		function iterateBackwards(obj, ishock)
			% iterate over the dynamic HJB, solver for policy functions, etc...

			% Parameters
			% ----------
			% ishock : the shock index, in reference to the shock sizes in the
			%	Params object used to instantiate this class
			%
			% Directly Modifies
			% -----------------
			% obj.cum_con_q1 : a cell array indexed by shock, where each entry
			%	contains the cumulative consumption over states given a shock in
			%	period 1 (where initial period is taken to be 0)
			%
			% obj.cum_con_q4 : a cell array indexed by shock, where each entry
			%	contains the cumulative consumption over states given a shock in
			%	period 4 (where initial period is taken to be 0)
			%
			% Note
			% ----
			% This method calls various other class methods and so indirectly
			% modifies other class properties.
			
            shock = obj.p.mpc_shocks(ishock);
			obj.cumcon = zeros(obj.p.nb_KFE*obj.p.na_KFE*obj.income.ny*obj.p.nz,4);

			dim = obj.p.nb_KFE*obj.p.na_KFE*obj.income.ny*obj.p.nz;

			for it = 4:-obj.p.delta_mpc:obj.p.delta_mpc
				timeUntilShock = obj.p.delta_mpc + 4 - it;
                timeUntilShock = round(timeUntilShock * 40) / 40;
				if mod(it,0.5) == 0
		            fprintf('\tUpdating cumulative consumption given news, quarter=%0.2f\n',it)
                end

                if obj.p.SDU == 1
			    	ez_adj = solver.SDU_income_risk_adjustment(obj.p, obj.V, obj.income);
			    end
                
                u_k = reshape(obj.KFEint.u,[],obj.p.ny);
                V_k = reshape(obj.V,[],obj.p.ny);
                V_k1 = zeros(size(V_k));
                for k = 1:obj.p.ny
                    ind1 = 1+obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*(k-1);
			    	ind2 = obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*k;
                    Ak = obj.A_HJB(ind1:ind2,ind1:ind2);
                    
                    indx_k = ~ismember(1:obj.income.ny,k);

                    if obj.p.SDU == 0
			    		Vk_stacked 	= sum(repmat(obj.income.ytrans(k,indx_k),obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz,1) ...
	                            .* V_k(:,indx_k),2);
			    		V_k1(:,k) = (obj.rho_diag + (1/obj.p.delta_mpc + obj.p.deathrate - obj.income.ytrans(k,k))...
		        				* speye(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz) - Ak)...
		                 	\ (u_k(:,k) + Vk_stacked + V_k(:,k)/obj.p.delta_mpc);
			    	else
			    		Vk_stacked 	= sum(squeeze(ez_adj(:, k, indx_k)) .* V_k(:,indx_k),2);
			    		V_k1(:,k) = (obj.rho_diag + (1/obj.p.delta_mpc + obj.p.deathrate)...
		        				* speye(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz) - Ak - ...
		        				spdiags(ez_adj(:, k, k), 0, obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz, obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz))...
		                 	\ (u_k(:,k) + Vk_stacked + V_k(:,k)/obj.p.delta_mpc);
			    	end
                        
                    
                end
                obj.V = reshape(V_k1,[obj.p.nb_KFE,obj.p.na_KFE,obj.income.ny,obj.p.nz]);

		        % find policies a fraction of a period back
		        obj.update_policies();
		        
		        % find A matrix a fraction of a period back
		        obj.update_A_matrix();

			    if obj.p.ComputeMPCS_news == 1
				    % matrix divisor for Feynman-Kac
				    FKmat = cell(1,obj.income.ny);
			        for k = 1:obj.income.ny
			        	ind1 = 1+obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*(k-1);
			        	ind2 = obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*k;

			        	if (obj.p.sigma_r > 0) && (obj.p.retrisk_KFE == 0)
			        		Ak = obj.A_FK(ind1:ind2, ind1:ind2);
			        	else
			        		Ak = obj.A_HJB(ind1:ind2, ind1:ind2);
			        	end

			        	FKmat{k} = speye(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz)*(...
				            	1/obj.p.delta_mpc + obj.p.deathrate - obj.income.ytrans(k,k)) - Ak;
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

		function update_policies(obj)
			% updates the policy functions during backward iteration

			% Modifies
			% --------
			% obj.KFEint : the intermediate policy functions over the
			%	KFE grids

			obj.KFEint = solver.find_policies(...
				obj.p, obj.income, obj.grids, obj.V);
		end

		function update_A_matrix(obj)
			% updates the transition matrix during backward iteration

			% Modifies
			% --------
			% obj.A_HJB : the transition matrix used for the policy functions,
			%	and when returns risk is ommitted from the KFE, is also used to
			%	solve the Feynman-Kac equations
			%
			% obj.A_FK : 
			obj.A_HJB = obj.A_constructor_HJB.construct(obj.KFEint, obj.V);

			if (obj.p.sigma_r > 0) && (obj.p.retrisk_KFE == 0)
				obj.A_FK = obj.A_constructor_FK.construct(obj.KFEint);
			end
		end

		function update_cum_con(obj, period, FKmat)
			% updates cumulative consumption for a subperiod during backward iteration

			% Parameters
			% ----------
			% period : the current quarter, integer
			%
			% FKmat : the matrix divisor for the Feynman-Kac equation
			%
			% Modifies
			% --------
			% obj.cumcon : the array of cumulative consumption, which is indexed
			%	by states and periods, of shape (nb_KFE*na_KFE*nz*ny, num_periods)

            cumcon_t_k = reshape(obj.cumcon(:,period),[],obj.p.ny);

			for k = 1:obj.income.ny
                ytrans_cc_k = sum(obj.ytrans_offdiag(k,:) .* cumcon_t_k,2);

                deathin_cc_k = obj.get_death_inflows(cumcon_t_k, k);

                ind1 = 1+obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*(k-1);
                ind2 = obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*k;
                RHS = reshape(obj.KFEint.c(:,:,:,k),[],1) + ytrans_cc_k + deathin_cc_k ...
                        + cumcon_t_k(:,k)/obj.p.delta_mpc;
                obj.cumcon(ind1:ind2,period) = FKmat{k}*RHS;
            end
		end

		function deathin_cc_k = get_death_inflows(obj, cumcon_t_k, k)
			% computes inflows for death

			% Parameters
			% ----------
			% cumcon_t_k : cumulative consumption for period t,
			%	of shape (nb_KFE*na_KFE*nz, ny)
			%
			% Returns
			% -------
			% deathin_cc_k : death inflows multiplied by cumulative consumption

            reshape_vec = [obj.p.nb_KFE*obj.p.na_KFE obj.p.nz obj.p.ny];
			cumcon_t_z_k = reshape(cumcon_t_k,reshape_vec);

            if (obj.p.Bequests == 1) && (obj.p.ResetIncomeUponDeath == 1)
                deathin_cc_k = obj.p.deathrate * sum(obj.income.ydist' .* cumcon_t_k(:,k),2);
                error('need to recode')
                if nz > 1
                    error('not correctly coded for nz > 1')
                end
            elseif (obj.p.Bequests == 1) && (obj.p.ResetIncomeUponDeath == 0)
                deathin_cc_k = obj.p.deathrate * cumcon_t_k(:,k);
            elseif (obj.p.Bequests == 0) && (obj.p.ResetIncomeUponDeath == 1)
                deathin_cc_k = obj.p.deathrate * sum(obj.income.ydist' .* cumcon_t_k(1,k),2);
                error('neet to recode')
                if nz > 1
                    error('not correctly coded for nz > 1')
                end
            elseif (obj.p.Bequests == 0) && (obj.p.ResetIncomeUponDeath == 0)
                deathin_cc_k = obj.p.deathrate * squeeze(cumcon_t_z_k(obj.grids.loc0b0a,:,k));
                deathin_cc_k = kron(deathin_cc_k,ones(obj.p.nb_KFE*obj.p.na_KFE,1));
            end
        end

		function computeMPCs(obj, pmf, ishock, cum_con_baseline)
			% performs MPC computations

            shock = obj.p.mpc_shocks(ishock);
            
            mpcs_1_quarterly = (obj.cum_con_q1{ishock} - cum_con_baseline(:,1)) / shock;
            obj.mpcs(ishock).avg_1_quarterly = mpcs_1_quarterly(:)' * pmf(:);

            mpcs_4_quarterly = (obj.cum_con_q4{ishock} - cum_con_baseline) / shock;
            obj.mpcs(ishock).avg_4_quarterly = mpcs_4_quarterly' * pmf(:);
            obj.mpcs(ishock).avg_4_annual = sum(obj.mpcs(ishock).avg_4_quarterly);
		end

		function savePolicies(obj,index,ishock)
			% saves the policy functions for MPC simulation
			
        	c = obj.KFEint.c;
        	s = obj.KFEint.s;
        	d = obj.KFEint.d;

	    	name = sprintf('policy%ishock%i.mat',index,ishock);
	    	save([obj.p.tempdirec name],'c','s','d')
        end
	end

end