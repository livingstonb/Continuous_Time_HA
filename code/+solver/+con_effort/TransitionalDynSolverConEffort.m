classdef TransitionalDynSolverConEffort < solver.TransitionalDynSolver
	% This class is used for solving for the policy functions
	% when a future shock is known, and optionally for computing
	% the MPCs out of news using Feynman-Kac.

	properties (SetAccess = protected)
		dim2Identity = 'c';
	end

	methods
		function obj = TransitionalDynSolverConEffort(params,income,grids,shocks)
			obj = obj@solver.TransitionalDynSolver(params,income,grids,shocks);
            obj.dim2 = params.nc_KFE;
		end

		function update_policies(obj)
			obj.KFEint.h = solver.con_effort.find_policies(obj.p,obj.income,obj.grids,obj.V);
            obj.KFEint.u = aux.u_fn(obj.KFEint.c,obj.p.riskaver) - aux.con_effort.con_adj_cost(obj.p,obj.KFEint.h)...
                - aux.con_effort.penalty(obj.grids.b.matrix,obj.p.penalty1,obj.p.penalty2);
            obj.KFEint.s = (obj.p.r_b+obj.p.deathrate*obj.p.perfectannuities) * obj.grids.b.matrix...
                + (1-obj.p.wagetax)* obj.income.y.matrixKFE - obj.KFEint.c;
        end

        function update_A_matrix(obj)
		    obj.A = solver.con_effort.construct_trans_matrix(obj.p,obj.income,obj.grids,obj.KFEint);
		end

		function deathin_cc_k = get_death_inflows(obj,cumcon_t,k)
			cumcon_t_k = reshape(cumcon_t,[],obj.income.ny);
            reshape_vec = [obj.p.nb_KFE*obj.p.nc_KFE*obj.p.nz obj.income.ny];
            cumcon_t_c_k = reshape(cumcon_t,reshape_vec);
            ytrans_cc_k = sum(obj.ytrans_offdiag(k,:) .* reshape(cumcon_t,[],obj.income.ny),2);

            if obj.p.Bequests == 1
                deathin_cc_k = obj.p.deathrate * cumcon_t_k(:,k);
            elseif obj.p.Bequests == 0
                deathin_cc_k = obj.p.deathrate * cumcon_t_c_k(1,:,k)';
                deathin_cc_k = kron(deathin_cc_k,ones(obj.p.nb_KFE,1));
            end
        end

        function savePolicies(obj,index,ishock)
        	h = obj.KFEint.h;
	    	name = sprintf('policy%ishock%i.mat',index,ishock);
	    	save([obj.p.tempdirec name],'h')
	    end
	end
end