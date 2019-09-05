classdef MPCSimulatorConEffort < statistics.MPCSimulator

	properties (SetAccess=protected)
		hinterp;
        hpolicy1;
        hpolicy2;
        
        timeIndex1;
        timeIndex2;
        
        savedTimesUntilShock;
	end

	methods
		function obj = MPCSimulatorConEffort(...
			p,income,grids,policies,shocks,shockperiod)
			obj = obj@statistics.MPCSimulator(...
				p,income,grids,policies,shocks,shockperiod,'c');

			if income.ny > 1
                obj.interp_grids = {grids.b.vec,grids.c.vec,income.y.vec};
            else
                obj.interp_grids = {grids.b.vec,grids.c.vec};
            end
            obj.hinterp{1} = griddedInterpolant(obj.interp_grids,policies.h,'makima');
            
            if shockperiod > 0
                % load times at which policy functions saved
                fname = [p.tempdirec 'savedTimesUntilShock.mat'];
                temp = load(fname);
                obj.savedTimesUntilShock = temp.savedTimesUntilShock;
            end
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
                if obj.shockperiod == 0
                    policyIndex = 1;
                else
                    policyIndex = ii;
                end
                
	            if income.ny > 1
	                hsim(:,ii) = ...
                        obj.hinterp{policyIndex}(obj.bsim(:,ii),obj.csim(:,ii),obj.ysim);
	            else
	                hsim(:,ii) = ...
	                	obj.hinterp{policyIndex}(obj.bsim(:,ii),obj.csim(:,ii));
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

		function update_interpolants(obj,p,time)
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
                    temp = load([p.tempdirec name]);
	    			obj.hpolicy1{ishock} = temp.h;
                    
	    			name = sprintf('policy%ishock%i.mat',obj.timeIndex2,ishock);
	    			temp = load([p.tempdirec name]);
	    			obj.hpolicy2{ishock} = temp.h;

	    			obj.hinterp{ii+1} = griddedInterpolant(...
	    				obj.interp_grids,obj.hpolicy1{ishock},'makima');
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
	    			temp = load([p.tempdirec name]);
                    obj.hpolicy1{ishock} = temp.h;
                    
                    if ~lastSavedTime
                        name = sprintf('policy%ishock%i.mat',obj.timeIndex2,ishock);
                        temp = load([p.tempdirec name]);
                        obj.hpolicy2{ishock} = temp.h;
                    end

	    			obj.hinterp{ii+1} = griddedInterpolant(...
	    				obj.interp_grids,obj.hpolicy1{ishock},'makima');
	    		end
	    	elseif timeUntilShock < obj.savedTimesUntilShock(end)
	    		timePastTime1 = time - obj.savedTimesUntilShock(obj.timeIndex1);
	    		timeDiff = obj.savedTimesUntilShock(obj.timeIndex2)...
	    			- obj.savedTimesUntilShock(obj.timeIndex1);
	    		proportion2 = timePastTime1 / timeDiff;

	    		for ii = 1:numel(obj.shocks)
                    ishock = obj.shocks(ii);
	    			h_approximate = (1-proportion2) * obj.hpolicy1{ishock}...
	    				+ proportion2 * obj.hpolicy2{ishock};
	    			obj.hinterp{ii+1} = griddedInterpolant(...
	    				obj.interp_grids,h_approximate,'makima');
	    		end
	    	elseif timeUntilShock > obj.savedTimesUntilShock(end)
	    		for ii = 1:numel(obj.shocks)
                    ishock = obj.shocks(ii);
	    			obj.hinterp{ii+1} = griddedInterpolant(...
	    				obj.interp_grids,obj.hpolicy1{ishock},'makima');
	    		end
	    	else
	    		error('logical error, hinterp failed to update')
	    	end
		end
	end
end