classdef MPCSimulatorConEffort < statistics.MPCSimulator

	properties (SetAccess=protected)
		hinterp;
	end

	methods
		function obj = MPCSimulatorConEffort(...
			p,income,grids,policies,shocks,shockperiod)
			obj = obj@statistics.MPCSimulator(...
				p,income,grids,policies,shocks,shockperiod,'c');

			if income.ny > 1
                interp_grids = {grids.b.vec,grids.c.vec,income.y.vec};
            else
                interp_grids = {grids.b.vec,grids.c.vec};
            end
            obj.hinterp{1} = griddedInterpolant(interp_grids,policies.h,'linear');
		end

		function draw_from_stationary_dist(obj,p,income,grids,pmf)

			index = draw_from_stationary_dist@statistics.MPCSimulator(...
				obj,p,income,grids,pmf);

			% initial consumption
			cgrid_flat = grids.c.matrix(:);
			obj.csim = repmat(cgrid_flat(index),1,obj.nshocks+1);
		end

	    function simulate_assets_one_period(obj,p,grids,current_csim)
	    	obj.bsim = obj.bsim + obj.mpc_delta ...
	            * ((p.r_b+p.deathrate*p.perfectannuities) * obj.bsim...
	                + (1-p.wagetax) * obj.ysim - current_csim);
	        obj.bsim = max(obj.bsim,grids.b.vec(1));
	        obj.bsim = min(obj.bsim,grids.b.vec(end));
	    end

	    function simulate_consumption_one_period(obj,p,income,grids)
	    	obj.cum_con = obj.cum_con + obj.csim * obj.mpc_delta;

	    	hsim = zeros(p.n_mpcsim,obj.nshocks+1);
	    	for ii = 1:obj.nshocks+1
	            if income.ny > 1
	                hsim(:,ii) = ...
	                	obj.hinterp{1}(obj.bsim(:,ii),obj.csim(:,ii),obj.ysim);
	            else
	                hsim(:,ii) = ...
	                	obj.hinterp{1}(obj.bsim(:,ii),obj.csim(:,ii));
	            end
	        end

	        if strcmp(p.hdef,'cdot')
	            obj.csim = obj.csim + obj.mpc_delta * hsim;
	        elseif strcmp(p.hdef,'cdot/c')
	            obj.csim = obj.csim + obj.mpc_delta * (hsim .* obj.csim);
	        end
	        obj.csim = max(obj.csim,grids.c.vec(1));
	        obj.csim = min(obj.csim,grids.c.vec(end));
	    end

		function update_interpolants(obj,p,time,shockperiod)
			if shockperiod == 0
				% policy function does not change in equilibrium
			else
				error('not yet coded')
			end
		end
	end
end