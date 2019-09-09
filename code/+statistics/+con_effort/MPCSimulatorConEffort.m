classdef MPCSimulatorConEffort < statistics.MPCSimulator
	% This class subclasses MPCSimulator to provide
	% MPC simulations for the consumption adj cost model.

	properties (SetAccess=protected)
		% interpolant for h = cdot or cdot/c policy
		hinterp;

		%  - the following are relevant only if simulating
		%  - MPCs out of news:

		% h policy function at beginning of current subperiod
        hpolicy1;

        % h policy function at end of current subperiod
        hpolicy2;
        
        % index of beginning of current supperiod within
        % saved policy functions
        timeIndex1;

        % index of end of current supperiod within
        % saved policy functions
        timeIndex2;
        
        % each element of this vector gives a time at
        % which a corresponding policy function was saved,
        % e.g. savedTimesUntilShock(5) is the time remaining
        % until the shock for the fifth saved policy function
        savedTimesUntilShock;
	end

	methods
		function obj = MPCSimulatorConEffort(...
			p,income,grids,policies,shocks,shockperiod)
			obj = obj@statistics.MPCSimulator(...
				p,income,grids,policies,shocks,shockperiod,'c');

			if (income.ny > 1) && (p.nz > 1)
                obj.interp_grids = {grids.b.vec,grids.c.vec,grids.z.vec,income.y.vec};
            elseif income.ny > 1
                obj.interp_grids = {grids.b.vec,grids.c.vec,income.y.vec};
            elseif p.nz > 1
                obj.interp_grids = {grids.b.vec,grids.c.vec,grids.z.vec};
            else
                obj.interp_grids = {grids.b.vec,grids.c.vec};
            end
            obj.hinterp{1} = griddedInterpolant(obj.interp_grids,squeeze(policies.h),'makima');
            
            if (shockperiod > 0) && (p.SimulateMPCS_news==1)
                % load times at which policy functions saved
                fname = [p.tempdirec 'savedTimesUntilShock.mat'];
                temp = load(fname);
                obj.savedTimesUntilShock = temp.savedTimesUntilShock;
            end

            for ishock = 1:6
            	obj.sim_mpcs(ishock).avg_quarterly_pos = NaN(4,1);
            	obj.sim_mpcs(ishock).avg_annual_pos = NaN;
                obj.sim_mpcs(ishock).responders_quarterly = NaN(4,1);
                obj.sim_mpcs(ishock).responders_annual = NaN;
            end
		end

		function draw_from_stationary_dist(obj,pmf)

			index = draw_from_stationary_dist@statistics.MPCSimulator(obj,pmf);

			% initial consumption
			cgrid_flat = obj.grids.c.matrix(:);
			obj.csim = repmat(cgrid_flat(index),1,obj.nshocks+1);
		end

	    function simulate_assets_one_period(obj,current_csim)
	    	obj.bsim = obj.bsim + obj.mpc_delta ...
	            * ((obj.p.r_b+obj.p.deathrate*obj.p.perfectannuities) * obj.bsim...
	                + (1-obj.p.wagetax) * obj.ysim - current_csim);
	        obj.bsim = max(obj.bsim,obj.grids.b.vec(1));
	        obj.bsim = min(obj.bsim,obj.grids.b.vec(end));
	    end

	    function simulate_consumption_one_period(obj)
	    	obj.cum_con = obj.cum_con + obj.csim * obj.mpc_delta;

	    	hsim = zeros(obj.p.n_mpcsim,obj.nshocks+1);
	    	for ii = 1:obj.nshocks+1
                if obj.shockperiod == 0
                    policyIndex = 1;
                else
                    policyIndex = ii;
                end
                
	            if obj.income.ny > 1
	                hsim(:,ii) = ...
                        obj.hinterp{policyIndex}(obj.bsim(:,ii),obj.csim(:,ii),obj.ysim);
	            else
	                hsim(:,ii) = ...
	                	obj.hinterp{policyIndex}(obj.bsim(:,ii),obj.csim(:,ii));
	            end
	        end

	        if strcmp(obj.p.hdef,'cdot')
	            obj.csim = obj.csim + obj.mpc_delta * hsim;
	        elseif strcmp(obj.p.hdef,'cdot/c')
	            obj.csim = obj.csim + obj.mpc_delta * (hsim .* obj.csim);
	        end
	        obj.csim = max(obj.csim,obj.grids.c.vec(1));
	        obj.csim = min(obj.csim,obj.grids.c.vec(end));
	    end

		function update_interpolants(obj,time)
           if obj.shockperiod == 0
               return
           end
            
            timeUntilShock = obj.shockperiod - time;
            
            if time == 0
	    		% load first two policy functions
	    		obj.timeIndex1 = find(obj.savedTimesUntilShock==obj.shockperiod);
	    		obj.timeIndex2 = obj.timeIndex1 + 1;

	    		for ii = 1:numel(obj.shocks)
                    ishock = obj.shocks(ii);
	    			name = sprintf('policy%ishock%i.mat',obj.timeIndex1,ishock);
                    temp = load([obj.p.tempdirec name]);
	    			obj.hpolicy1{ishock} = temp.h;
                    
	    			name = sprintf('policy%ishock%i.mat',obj.timeIndex2,ishock);
	    			temp = load([obj.p.tempdirec name]);
	    			obj.hpolicy2{ishock} = temp.h;

	    			obj.hinterp{ii+1} = griddedInterpolant(...
	    				obj.interp_grids,squeeze(obj.hpolicy1{ishock}),'makima');
	    		end
	    	elseif abs(timeUntilShock - obj.savedTimesUntilShock(obj.timeIndex2)) < 1e-6
	    		% update policy functions
	    		obj.timeIndex1 = obj.timeIndex1 + 1;
                
                if obj.timeIndex2 == numel(obj.savedTimesUntilShock)
                    lastSavedTime = true;
                else
                    lastSavedTime = false;
                    obj.timeIndex2 = obj.timeIndex2 + 1;
                end

	    		for ii = 1:numel(obj.shocks)
                    ishock = obj.shocks(ii);
	    			name = sprintf('policy%ishock%i.mat',obj.timeIndex1,ishock);
	    			temp = load([obj.p.tempdirec name]);
                    obj.hpolicy1{ishock} = temp.h;
                    
                    if ~lastSavedTime
                        name = sprintf('policy%ishock%i.mat',obj.timeIndex2,ishock);
                        temp = load([obj.p.tempdirec name]);
                        obj.hpolicy2{ishock} = temp.h;
                    end

	    			obj.hinterp{ii+1} = griddedInterpolant(...
	    				obj.interp_grids,squeeze(obj.hpolicy1{ishock}),'makima');
	    		end
	    	elseif timeUntilShock > obj.savedTimesUntilShock(obj.timeIndex2)
	    		timeUntil2 = timeUntilShock - obj.savedTimesUntilShock(obj.timeIndex2);
	    		timeDiff = obj.savedTimesUntilShock(obj.timeIndex1)...
	    			- obj.savedTimesUntilShock(obj.timeIndex2);
	    		proportion1 = timeUntil2 / timeDiff;

	    		for ii = 1:numel(obj.shocks)
                    ishock = obj.shocks(ii);
	    			h_approximate = proportion1 * obj.hpolicy1{ishock}...
	    				+ (1-proportion1) * obj.hpolicy2{ishock};
	    			obj.hinterp{ii+1} = griddedInterpolant(...
	    				obj.interp_grids,squeeze(h_approximate),'makima');
	    		end
	    	elseif timeUntilShock < obj.savedTimesUntilShock(end)
	    		for ii = 1:numel(obj.shocks)
                    ishock = obj.shocks(ii);
	    			obj.hinterp{ii+1} = griddedInterpolant(...
	    				obj.interp_grids,squeeze(obj.hpolicy1{ishock}),'makima');
	    		end
	    	else
	    		error('logical error, hinterp failed to update')
	    	end
		end

		function compute_mpcs(obj,period)
			compute_mpcs@statistics.MPCSimulator(obj,period);

		    for is = 1:obj.nshocks
		    	ishock = obj.shocks(is); % index of shock in parameters
		    	shock = obj.p.mpc_shocks(ishock);

		    	con_diff = obj.shock_cum_con{ishock} - obj.baseline_cum_con;
            
            	if (shock > 0) || ismember(obj.shockperiod,[0,1])
	            	obj.sim_mpcs(ishock).avg_quarterly_pos(period) = mean(con_diff(con_diff(:,period)/shock>0,period) / shock);
					obj.sim_mpcs(ishock).responders_quarterly(period) = mean(con_diff(:,period)/shock>0);
	            end

	            if (shock > 0) || (obj.shockperiod == 4)
                    obj.sim_mpcs(ishock).avg_annual_pos = mean(sum(con_diff(con_diff/shock>0),2) / shock);
                    obj.sim_mpcs(ishock).responders_annual = mean(sum(con_diff/shock,2)>0);
                end
		    end
		end
	end
end