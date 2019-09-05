classdef MPCSimulatorTwoAsset < statistics.MPCSimulator
	% This class subclasses MPCSimulator to provide
	% MPC simulations for the two-asset model.

	properties (SetAccess=protected)
		% illiquid asset
        asim;

        % policy function interpolants
		dinterp;
		sinterp;
		cinterp;
	end

	methods
		function obj = MPCSimulatorTwoAsset(...
			p,income,grids,policies,shocks,nperiods)
			obj = obj@statistics.MPCSimulator(...
				p,income,grids,policies,shocks,nperiods,'a');

			if income.ny > 1
                interp_grids = {obj.grids.b.vec,obj.grids.a.vec,obj.income.y.vec};
            else
                interp_grids = {obj.grids.b.vec,obj.grids.a.vec};
            end
            obj.dinterp = griddedInterpolant(interp_grids,policies.d,'linear');
			obj.cinterp = griddedInterpolant(interp_grids,policies.c,'linear');
			obj.sinterp = griddedInterpolant(interp_grids,policies.s,'linear');
		end

		function draw_from_stationary_dist(obj,pmf)
			% this function draws from the stationary distribution
			index = draw_from_stationary_dist@statistics.MPCSimulator(obj,pmf);

			% initial illiquid assets
			agrid_flat = obj.grids.a.matrix(:);
			obj.asim = repmat(agrid_flat(index),1,obj.nshocks+1);
		end

	    function simulate_assets_one_period(obj,~)
	    	% this function simulates assets over the next time
	    	% delta

	    	% interpolate to find decisions
	    	s = obj.sinterp(obj.bsim(:),obj.asim(:),obj.yrep);
	    	d = obj.dinterp(obj.bsim(:),obj.asim(:),obj.yrep);

	    	s = reshape(s,[],obj.nshocks+1);
	    	d = reshape(d,[],obj.nshocks+1);

	    	% update liquid assets
	    	obj.bsim = obj.bsim + obj.mpc_delta ...
				* (s - d - aux.two_asset.adj_cost_fn(d,obj.asim,obj.p));
			obj.bsim = max(obj.bsim,obj.grids.b.vec(1));
			obj.bsim = min(obj.bsim,obj.grids.b.vec(end));

			% update illiquid assets
	    	obj.asim = obj.asim + obj.mpc_delta ...
	            * (d + (obj.p.r_a+obj.p.deathrate*obj.p.perfectannuities) * obj.asim ...
				+ obj.p.directdeposit*obj.ysim);
	        obj.asim = max(obj.asim,obj.grids.a.vec(1));
	        obj.asim = min(obj.asim,obj.grids.a.vec(end));
	    end

	    function simulate_consumption_one_period(obj)
	    	c = obj.cinterp(obj.bsim(:),obj.asim(:),obj.yrep);
	    	c = reshape(c,[],obj.nshocks+1);

	    	obj.cum_con = obj.cum_con + c * obj.mpc_delta;
	    end
	end
end