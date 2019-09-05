classdef MPCSimulatorTwoAsset < statistics.MPCSimulator

	properties (SetAccess=protected)
        asim;

		dinterp;
		sinterp;
		cinterp;
	end

	methods
		function obj = MPCSimulatorTwoAsset(...
			p,income,grids,policies,shocks,nperiods,dim2Identity)
			obj = obj@statistics.MPCSimulator(...
				p,income,grids,policies,shocks,nperiods,dim2Identity);

			if income.ny > 1
                interp_grids = {grids.b.vec,grids.a.vec,income.y.vec};
            else
                interp_grids = {grids.b.vec,grids.a.vec};
            end
            obj.dinterp = griddedInterpolant(interp_grids,policies.d,'linear');
			obj.cinterp = griddedInterpolant(interp_grids,policies.c,'linear');
			obj.sinterp = griddedInterpolant(interp_grids,policies.s,'linear');
		end

		function draw_from_stationary_dist(obj,p,income,grids,pmf)

			index = draw_from_stationary_dist@statistics.MPCSimulator(...
				obj,p,income,grids,pmf);

			% initial illiquid assets
			agrid_flat = grids.a.matrix(:);
			obj.asim = repmat(agrid_flat(index),1,obj.nshocks+1);
		end

	    function simulate_assets_one_period(obj,p,grids,~)
	    	s = obj.sinterp(obj.bsim(:),obj.asim(:),obj.yrep);
	    	d = obj.dinterp(obj.bsim(:),obj.asim(:),obj.yrep);

	    	s = reshape(s,[],obj.nshocks+1);
	    	d = reshape(d,[],obj.nshocks+1);

	    	obj.bsim = obj.bsim + obj.mpc_delta ...
				* (s - d - aux.two_asset.adj_cost_fn(d,obj.asim,p));
			obj.bsim = max(obj.bsim,grids.b.vec(1));
			obj.bsim = min(obj.bsim,grids.b.vec(end));

	    	obj.asim = obj.asim + obj.mpc_delta ...
	            * (d + (p.r_a+p.deathrate*p.perfectannuities) * obj.asim ...
				+ p.directdeposit*obj.ysim);
	        obj.asim = max(obj.asim,grids.a.vec(1));
	        obj.asim = min(obj.asim,grids.a.vec(end));
	    end

	    function simulate_consumption_one_period(obj,~,~,~)
	    	c = obj.cinterp(obj.bsim(:),obj.asim(:),obj.yrep);
	    	c = reshape(c,[],obj.nshocks+1);

	    	obj.cum_con = obj.cum_con + c * obj.mpc_delta;
	    end
	end
end