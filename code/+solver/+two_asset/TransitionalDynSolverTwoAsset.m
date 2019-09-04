classdef TransitionalDynSolverTwoAsset < solver.TransitionalDynSolver

	properties (SetAccess = protected)
		dim2Identity = 'a';
	end

	methods
		function obj = TransitionalDynSolverTwoAsset(params,income,grids)
			obj = obj@solver.TransitionalDynSolver(params,income,grids);
            obj.dim2 = params.na_KFE;
		end

		function update_policies(obj)
            V_normaldim = permute(obj.V,[1 2 4 3]);
			obj.KFEint = solver.two_asset.find_policies(...
				obj.p,obj.income,obj.grids,V_normaldim);
		end

		function update_A_matrix(obj)
		    obj.A = solver.two_asset.construct_trans_matrix(...
		    	obj.p,obj.income,obj.grids,obj.KFEint,'KFE');
		end

		function deathin_cc_k = get_death_inflows(obj,cumcon_t,k)
			cumcon_t_k = reshape(cumcon_t,[],obj.p.ny);
            reshape_vec = [obj.p.nb_KFE*obj.p.na_KFE obj.p.nz obj.p.ny];
			cumcon_t_z_k = reshape(cumcon_t_k,reshape_vec);
            ytrans_cc_k = sum(obj.ytrans_offdiag(k,:) .* reshape(cumcon_t,[],obj.p.ny),2);

            if (obj.p.Bequests == 1) && (obj.p.ResetIncomeUponDeath == 1)
                deathin_cc_k = obj.p.deathrate * sum(obj.income.ydist' .* cumcon_t_k,2);
                if nz > 1
                    error('not correctly coded for nz > 1')
                end
            elseif (obj.p.Bequests == 1) && (obj.p.ResetIncomeUponDeath == 0)
                deathin_cc_k = obj.p.deathrate * cumcon_t_k(:,k);
            elseif (obj.p.Bequests == 0) && (obj.p.ResetIncomeUponDeath == 1)
                deathin_cc_k = obj.p.deathrate * sum(obj.income.ydist' .* cumcon_t_k(1,:),2);
                if nz > 1
                    error('not correctly coded for nz > 1')
                end
            elseif (obj.p.Bequests == 0) && (obj.p.ResetIncomeUponDeath == 0)
                deathin_cc_k = obj.p.deathrate * cumcon_t_z_k(1,:,k)';
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