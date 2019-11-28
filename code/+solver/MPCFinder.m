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
		% p : Model parameters
		p; % parameters
		income;
		grids;
        
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
		function obj = MPCFinder(p, income, grids)
			% class constructor

			% Parameters
			% ----------
			% p : a Params object
			% 
			% income : an Income object
			%
			% grids : a Grid object, providing the asset grids
			%	for the KFE
			%
			% Returns
			% -------
			% obj : an MPCFinder object

			obj.p = p;
			obj.income = income;
			obj.grids = grids;

			obj.ResetIncomeUponDeath = p.ResetIncomeUponDeath;

			obj.ytrans_offdiag = income.ytrans - diag(diag(income.ytrans));

			for ii = 1:6
				obj.mpcs(ii).mpcs = NaN;
				obj.mpcs(ii).quarterly = NaN(4,1);
                obj.mpcs(ii).annual = NaN;
			end
		end

		function solve(obj, KFE, pmf, A)
			% computes the MPCs using Feynman-Kac by calling
			% a sequence of class methods

			% Parameters
			% ----------
			% KFE : a structure containing the policy functions
			%	on the KFE grids
			%
			% pmf : the equilibrium probability mass function,
			%	of shape (nb_KFE, na_KFE, nz, ny)
			%
			% A : the transition matrix for the KFE
			%
			% Directly Modifies
			% -----------------
			% obj.mpcs : the computed MPC statistics
			%
			% Note
			% ----
			% This method also modifies other class properties
			% via other methods.

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
	end

	methods (Access=private)
		function create_FK_matrices(obj, A)
			% finds the matrix divisor for the Feynman-Kac formula

			% Parameters
			% ----------
			% A : the KFE transition matrix
			%
			% Modifies
			% --------
			% obj.FKmat: a cell array containing the FK matrix divisors
			%	for each income block, shape (1, ny)

			obj.FKmat = cell(1,obj.income.ny);
		    for k = 1:obj.income.ny
		    	ind1 = 1+obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*(k-1);
		    	ind2 = obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*k;

		        obj.FKmat{k} = speye(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz)*(1/obj.p.delta_mpc ...
		        					+ obj.p.deathrate - obj.income.ytrans(k,k))...
		                		- A(ind1:ind2,ind1:ind2);
		        obj.FKmat{k} = inverse(obj.FKmat{k});
		    end
		end

		function iterate_backward(obj, KFE)
			% iterates backward four quarters on the Feynman-Kac equation
			% to compute cum_con_baseline
			%
			% each quarter is split into 1/delta_mpc subperiods, which
			% are then iterated over

			% Parameters
			% ----------
			% KFE : a structure containing the policy functions and the
			%	value function, on the KFE grid
			%
			% Directly Modifies
			% -----------------
			% obj.cum_con_baseline : cumulative consumption at a given
			%	quarter for the baseline case, where the second dimension
			%	is the quarter; of shape (nb_KFE*na_KFE*nz*ny, num_quarters)
			%
			% Note
			% ----
			% This method also modifies other class properties
			% via other methods.

			dim = obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*obj.income.ny;
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

			reshape_vec = [obj.p.nb_KFE,obj.p.na_KFE,obj.p.nz,obj.income.ny,4];
		    obj.cum_con_baseline = reshape(obj.cum_con_baseline,[],4);
		end

		function update_cumcon(obj, KFE, period)
			% update cumulative consumption in the baseline case by
			% iterating on the Feynman-Kac equation

			% Parameters
			% ----------
			% KFE : a structure containing the policy functions and the
			%	value function, on the KFE grid
			%
			% period : the current period, integer
			%
			% Modifies
			% --------
			% obj.cumcon : cumulative consumption, modified for the given
			%	period, and has shape (nb_KFE*na_KFE*nz*ny, num_periods)

			cumcon_t = obj.cumcon(:,period);
			cumcon_t_k = reshape(cumcon_t,[],obj.income.ny);

            reshape_vec = [obj.p.nb_KFE*obj.p.na_KFE obj.p.nz obj.income.ny];
            cumcon_t_z_k = reshape(cumcon_t_k,reshape_vec);
                
			for k = 1:obj.income.ny
				ytrans_cc_k = sum(obj.ytrans_offdiag(k,:) .* cumcon_t_k,2);
  
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
                    deathin_cc_k = obj.p.deathrate * cumcon_t_z_k(obj.grids.loc0b0a,:,k)';
                    deathin_cc_k = kron(deathin_cc_k,ones(obj.p.nb_KFE*obj.p.na_KFE,1));
                end

                ind1 = 1+obj.p.na_KFE*obj.p.nb_KFE*obj.p.nz*(k-1);
                ind2 = obj.p.na_KFE*obj.p.nb_KFE*obj.p.nz*k;
                RHS = reshape(KFE.c(:,:,:,k),[],1) + ytrans_cc_k + deathin_cc_k ...
                            + cumcon_t_k(:,k)/obj.p.delta_mpc;
                obj.cumcon(ind1:ind2,period) = obj.FKmat{k}*RHS;
			end
		end

		function cumulative_consumption_with_shock(obj, ishock)
			% use baseline cumulative consumption to approximate
			% cumulative consumption for households presented with
			% an income shock in the first period

			% Parameters
			% ----------
			% ishock : the index of the shock, in reference to the
			%	shocks vector contained in the Params object used to
			%	instantiate this class
			%
			% Modifies
			% --------
			% obj.cum_con_shock : a cell array, indexed by shock, containing
			%	cumulative consumption over states for a given period; each
			%	cell contains an array of shape (nb_KFE*na_KFE*nz*ny, num_periods)

			shock = obj.p.mpc_shocks(ishock);
			bgrid_mpc_vec = obj.grids.b.vec + shock;

			if shock < 0
	            below_bgrid = bgrid_mpc_vec < obj.grids.b.vec(1);
	            bgrid_mpc_vec(below_bgrid) = obj.grids.b.vec(1);
	        end

	        dim = obj.p.na_KFE*obj.income.ny*obj.p.nz;
	        bgrid_mpc = repmat(bgrid_mpc_vec,dim,1);

	        % grids for interpolation
	        interp_grids = {obj.grids.b.vec};
            interp_grids{end+1} = obj.grids.a.vec;


	        if obj.p.nz > 1
	        	interp_grids{end+1} = obj.grids.z.vec;
	        end

	        if obj.income.ny > 1
	        	interp_grids{end+1} = obj.income.y.vec;
	        end

			reshape_vec = [obj.p.nb_KFE obj.p.na_KFE obj.p.nz obj.income.ny];
			for period = 1:4
				% cumulative consumption in 'period'
	            con_period = reshape(obj.cum_con_baseline(:,period),reshape_vec);
	            mpcinterp = griddedInterpolant(interp_grids,squeeze(con_period),'linear');

	            if (obj.income.ny > 1) && (obj.p.nz > 1)
	                obj.cum_con_shock{ishock}(:,period) = mpcinterp(bgrid_mpc,obj.grids.a.matrix(:),obj.grids.z.matrix(:),obj.income.y.matrixKFE(:));
	            elseif (obj.income.ny==1) && (obj.p.nz > 1)
	                obj.cum_con_shock{ishock}(:,period) = mpcinterp(bgrid_mpc,obj.grids.a.matrix(:),obj.grids.z.matrix(:));
	            elseif obj.income.ny > 1
	                obj.cum_con_shock{ishock}(:,period) = mpcinterp(bgrid_mpc,obj.grids.a.matrix(:),obj.income.y.matrixKFE(:));
	            else
	                obj.cum_con_shock{ishock}(:,period) = mpcinterp(bgrid_mpc,obj.grids.a.matrix(:));
	            end

	            if (shock < 0) && (sum(below_bgrid)>0) && (period==1)
	                temp = reshape(obj.cum_con_shock{ishock}(:,period),reshape_vec);
	                temp(below_bgrid,:,:,:) = con_period(1,:,:,:) + shock + obj.grids.b.vec(below_bgrid) - obj.grids.b.vec(1);
	                obj.cum_con_shock{ishock}(:,period) = temp(:);                      
	            end
	        end
		end

		function computeMPCs(obj, pmf, ishock)
			% compute MPCs using the cumulative consumption arrays
			% found previously

			% Parameters
			% ----------
			% pmf : the equilibrium probability mass function of the
			%	baseline, of shape (nb_KFE, na_KFE, nz, ny)
			%
			% ishock : the shock index, in reference to the shock vector
			%	in the Params object used to instantiate this class
			%
			% Modifies
			% --------
			% obj.mpcs : the final MPC statistics computed from this class,
			%	a structure array of size nshocks
			
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