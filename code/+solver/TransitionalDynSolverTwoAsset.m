classdef TransitionalDynSolverTwoAsset < solver.TransitionalDynSolver
	% This class is used for solving for the policy functions
	% when a future shock is known, and optionally for computing
	% the MPCs out of news using Feynman-Kac.

	properties (SetAccess = protected)
		dim2Identity = 'a';
	end

	methods
		function obj = TransitionalDynSolverTwoAsset(params,income,grids,shocks)
			obj = obj@solver.TransitionalDynSolver(params,income,grids,shocks);
            obj.dim2 = params.na_KFE;

            blockDim = params.nb_KFE*params.na_KFE*params.nz;
			if numel(params.rhos) > 1
		        rhocol = kron(params.rhos',ones(params.nb_KFE*params.na_KFE,1));
		        obj.rho_diag = spdiags(rhocol,0,blockDim,blockDim);
		    else
		        obj.rho_diag = params.rho * speye(blockDim);
			end
		end

		function update_policies(obj)
			obj.KFEint = solver.find_policies(...
				obj.p,obj.income,obj.grids,obj.V);
		end

		function update_A_matrix(obj)
			obj.A_HJB = obj.A_constructor_HJB.construct(obj.KFEint, obj.V);

			if (obj.p.sigma_r > 0) && (obj.p.retrisk_KFE == 0)
				obj.A_FK = obj.A_constructor_FK.construct(obj.KFEint);
			end
		end

		function deathin_cc_k = get_death_inflows(obj,cumcon_t_k,k)
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

        function savePolicies(obj,index,ishock)
        	c = obj.KFEint.c;
        	s = obj.KFEint.s;
        	d = obj.KFEint.d;

	    	name = sprintf('policy%ishock%i.mat',index,ishock);
	    	save([obj.p.tempdirec name],'c','s','d')
        end
	end

end