classdef MPCFinder < handle
	% This class contains methods that compute MPCs
	% over certain time periods. Baseline expected
	% consumption is found, and then this quantity--via
	% interpolation--is used to find expected consumption 
	% given an asset shock. Then MPCs are computed.
	%
	% Once the class is instantiated, the 'mpcs' results structure
	% is populated with NaN's. Calling solve() will compute
	% the MPCs.

	properties (SetAccess=private)
		p; % parameters
		income;
		grids;
        
		dim2;
		dim2Identity;
        
        ResetIncomeUponDeath;

		FKmat; % Feynman-Kac divisor matrix

		% cumulative consumption
		cumcon; % current state
		cum_con_baseline; % baseline
		cum_con_shock = cell(1,6); % shocked

		% income transitions w/o diagonal
		ytrans_offdiag;

		% results structure
		mpcs = struct();

		solved = false;
	end

	methods
		function obj = MPCFinder(p,income,grids,dim2Identity)
			obj.p = p;
			obj.income = income;
			obj.grids = grids;

			obj.dim2Identity = dim2Identity;

			if strcmp(dim2Identity,'a')
				obj.dim2 = p.na_KFE;
				obj.ResetIncomeUponDeath = p.ResetIncomeUponDeath;
			elseif strcmp(dim2Identity,'c')
				obj.dim2 = p.nc_KFE;
				obj.ResetIncomeUponDeath = 0;
			end

			obj.ytrans_offdiag = income.ytrans - diag(diag(income.ytrans));

			for ii = 1:6
				obj.mpcs(ii).mpcs = NaN;
				obj.mpcs(ii).quarterly = NaN(4,1);
			end
		end

		function solve(obj,KFE,pmf,A)
			% This is the function to call after instantiating the class,
			% to solve for the MPCs. The variable 'pmf' is the pmf over
			% the (b,a,y,z)-space on the KFE grids and Au is the transition
			% matrix on the KFE grids.

			if obj.solved
				error('Already solved, create another instance')
			end

			obj.create_FK_matrices(A);
			obj.iterate_backward(KFE);

			for ishock = 1:6
				obj.cumulative_consumption_with_shock(ishock);
				obj.computeMPCs(pmf,ishock)
			end

			obj.solved = true;
		end

		function create_FK_matrices(obj,A)
			% finds the matrix divisor in the Feynman-Kac formula
			obj.FKmat = cell(1,obj.income.ny);
		    for k = 1:obj.income.ny
		    	ind1 = 1+obj.p.nb_KFE*obj.dim2*obj.p.nz*(k-1);
		    	ind2 = obj.p.nb_KFE*obj.dim2*obj.p.nz*k;
		        obj.FKmat{k} = speye(obj.p.nb_KFE*obj.dim2*obj.p.nz)*(1/obj.p.delta_mpc ...
		        					+ obj.p.deathrate - obj.income.ytrans(k,k))...
		                		- A(ind1:ind2,ind1:ind2);
		        obj.FKmat{k} = inverse(obj.FKmat{k});
		    end
		end

		function iterate_backward(obj,KFE)
			% iterates backward four quarters on the Feynman-Kac equation
			% to compute cum_con_baseline
			%
			% each quarter is split into 1/delta_mpc subperiods, which
			% are then iterated over

			dim = obj.p.nb_KFE*obj.dim2*obj.p.nz*obj.income.ny;
			obj.cumcon = zeros(dim,4);
			for it = 4:-obj.p.delta_mpc:obj.p.delta_mpc
				if mod(it*4,1) == 0
					fprintf('\tUpdating baseline cumulative consumption, quarter=%0.2f\n',it)
				end

				for period = ceil(it):4
					% when 'it' falls to 'period', start updating
					% that 'period'
					obj.update_cumcon(KFE,period)
				end
			end

			obj.cum_con_baseline = zeros(dim,4);
			obj.cum_con_baseline(:,1) = obj.cumcon(:,1);
			for period = 2:4
				obj.cum_con_baseline(:,period) = obj.cumcon(:,period)...
					- obj.cumcon(:,period-1);
			end

			reshape_vec = [obj.p.nb_KFE,obj.dim2,obj.p.nz,obj.income.ny,4];
			obj.cum_con_baseline = reshape(obj.cum_con_baseline,reshape_vec);
		    obj.cum_con_baseline = permute(obj.cum_con_baseline,[1 2 4 3 5]);
		    obj.cum_con_baseline = reshape(obj.cum_con_baseline,[],4);
		end

		function update_cumcon(obj,KFE,period)
			% update cumulative consumption using the Feynman-Kac
			% iterative procedure
			for k = 1:obj.income.ny
				cumcon_t = obj.cumcon(:,period);
				cumcon_t_k = reshape(cumcon_t,[],obj.income.ny);

				if strcmp(obj.dim2Identity,'a')
					reshape_vec = [obj.p.nb_KFE*obj.dim2 obj.p.nz obj.income.ny];
					cumcon_t_z_k = reshape(cumcon_t_k,reshape_vec);
				elseif strcmp(obj.dim2Identity,'c')
					reshape_vec = [obj.p.nb_KFE obj.dim2 obj.income.ny];
					cumcon_t_c_k = reshape(cumcon_t_k,reshape_vec);
				end

                ytrans_cc_k = sum(obj.ytrans_offdiag(k,:) .* reshape(cumcon_t,[],obj.income.ny),2);

                if (obj.p.Bequests == 1) && (obj.ResetIncomeUponDeath == 1)
                	deathin_cc_k = obj.p.deathrate * sum(obj.income.ydist' .* cumcon_t_k,2);
                    if nz > 1
                        error('not correctly coded for nz > 1')
                    end
                elseif (obj.p.Bequests == 1) && (obj.ResetIncomeUponDeath == 0)
                	deathin_cc_k = obj.p.deathrate * cumcon_t_k(:,k);
                elseif (obj.p.Bequests == 0) && (obj.ResetIncomeUponDeath == 1)
                	deathin_cc_k = obj.p.deathrate * sum(obj.income.ydist' .* cumcon_t_k(1,:),2);
                    if nz > 1
                        error('not correctly coded for nz > 1')
                    end
                elseif (obj.p.Bequests == 0) && (obj.ResetIncomeUponDeath == 0)
                	if strcmp(obj.dim2Identity,'a')
	                	deathin_cc_k = obj.p.deathrate * cumcon_t_z_k(1,:,k)';
	                    deathin_cc_k = kron(deathin_cc_k,ones(obj.p.nb_KFE*obj.dim2,1));
	                elseif strcmp(obj.dim2Identity,'c')
	                    deathin_cc_k = obj.p.deathrate * cumcon_t_c_k(1,:,k);
	                    deathin_cc_k = kron(deathin_cc_k,ones(obj.dim2,1));
	                end
                end

                ind1 = 1+obj.dim2*obj.p.nb_KFE*obj.p.nz*(k-1);
                ind2 = obj.dim2*obj.p.nb_KFE*obj.p.nz*k;
                RHS = reshape(KFE.c(:,:,k,:),[],1) + ytrans_cc_k + deathin_cc_k ...
                            + cumcon_t(ind1:ind2)/obj.p.delta_mpc;
                obj.cumcon(ind1:ind2,period) = obj.FKmat{k}*RHS;
			end
		end

		function cumulative_consumption_with_shock(obj,ishock)
			% interpolate cumulative consumption with a shock to
			% assets
			shock = obj.p.mpc_shocks(ishock);
			bgrid_mpc_vec = obj.grids.b.vec + shock;

			if shock < 0
	            below_bgrid = bgrid_mpc_vec < obj.grids.b.vec(1);
	            bgrid_mpc_vec(below_bgrid) = obj.grids.b.vec(1);
	        end

	        dim = obj.dim2*obj.income.ny*obj.p.nz;
	        bgrid_mpc = repmat(bgrid_mpc_vec,dim,1);

	        % grids for interpolation
	        interp_obj.grids = {obj.grids.b.vec};
	        if strcmp(obj.dim2Identity,'a')
	        	interp_obj.grids{end+1} = obj.grids.a.vec;
	        elseif strcmp(obj.dim2Identity,'c')
	        	interp_obj.grids{end+1} = obj.grids.c.vec;
	        end

	        if obj.income.ny > 1
	        	interp_obj.grids{end+1} = obj.income.y.vec;
	        end

	        if obj.p.nz > 1
	        	interp_obj.grids{end+1} = obj.p.rhos';
	        end

	        if strcmp(obj.dim2Identity,'a')
	        	dim2_mat = obj.grids.a.matrix(:);
	        	if numel(obj.p.rhos) > 1
		        	dim = obj.p.nb_KFE*obj.dim2*obj.income.ny;
		        	heterog_mpcs = kron(obj.p.rhos',ones(dim,1));
				end
			else
				dim2_mat = obj.grids.c.matrix(:);
			end

			reshape_vec = [obj.p.nb_KFE,obj.dim2,obj.income.ny,obj.p.nz];
			for period = 1:4
				% cumulative consumption in 'period'
	            con_period = reshape(obj.cum_con_baseline(:,period),reshape_vec);
	            con_period = squeeze(con_period);
	            mpcinterp = griddedInterpolant(interp_obj.grids,con_period,'linear');
	            con_period = reshape(con_period,reshape_vec);

	            if (obj.income.ny > 1) && (obj.p.nz > 1)
	                obj.cum_con_shock{ishock}(:,period) = mpcinterp(bgrid_mpc,dim2_mat,obj.income.y.matrixKFE(:),heterog_mpcs);
	            elseif obj.p.nz > 1
	                obj.cum_con_shock{ishock}(:,period) = mpcinterp(bgrid_mpc,dim2_mat,heterog_mpcs);
	            elseif obj.income.ny > 1
	                obj.cum_con_shock{ishock}(:,period) = mpcinterp(bgrid_mpc,dim2_mat,obj.income.y.matrixKFE(:));
	            else
	                obj.cum_con_shock{ishock}(:,period) = mpcinterp(bgrid_mpc,dim2_mat);
	            end

	            if (shock < 0) && (sum(below_bgrid)>0) && (period==1)
	                temp = reshape(obj.cum_con_shock{ishock}(:,period),reshape_vec);
	                temp(below_bgrid,:,:,:) = con_period(1,:,:,:) + shock + obj.grids.b.vec(below_bgrid);
	                obj.cum_con_shock{ishock}(:,period) = temp(:);                      
	            end
	        end
		end

		function computeMPCs(obj,pmf,ishock)
			% compute MPCs from cumulative consumption
			shock = obj.p.mpc_shocks(ishock);

			% MPCs out of a shock at beginning of quarter 0
			mpcs = (obj.cum_con_shock{ishock} - obj.cum_con_baseline) / shock;
			if ishock == 5
				obj.mpcs(5).mpcs = mpcs;
			end

			obj.mpcs(ishock).quarterly = mpcs' * pmf(:);
			obj.mpcs(ishock).annual = sum(mpcs,2)' * pmf(:);
		end
	end
end