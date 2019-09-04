classdef RhoBoundsFinder < handle

	properties (SetAccess = private)
		neg_stepsize;
		pos_stepsize;
		rho_bound_history = [];
		current_rho_bound;
		lastDirection = 'none';
		lag;
	end

	methods
		function obj = RhoBoundsFinder(neg_stepsize,pos_stepsize,rho0)
			obj.neg_stepsize = neg_stepsize;
			obj.pos_stepsize = pos_stepsize;
			obj.current_rho_bound = rho0;
			obj.rho_bound_history(end+1) = obj.current_rho_bound;
		end

        function increase(obj)
        	disp('Raising candidate for rho bound')
        	if strcmp(obj.lastDirection,'none') || strcmp(obj.lastDirection,'up')
        		% increase by step
        		obj.current_rho_bound = obj.current_rho_bound + obj.pos_stepsize;
        		obj.lag = 1;
        	elseif strcmp(obj.lastDirection,'down')
        		% update by midpoint formula
        		obj.current_rho_bound = (obj.rho_bound_history(end-obj.lag)...
        			+ obj.current_rho_bound)/2;
        		obj.lag = obj.lag + 1;
        	end

        	obj.rho_bound_history(end+1) = obj.current_rho_bound;
        	obj.lastDirection = 'up';
        end

        function decrease(obj)
        	disp('Reducing candidate for rho bound')
        	if strcmp(obj.lastDirection,'none') || strcmp(obj.lastDirection,'down')
        		% increase by step
        		obj.current_rho_bound = obj.current_rho_bound + obj.neg_stepsize;
        		obj.lag = 1;
        	elseif strcmp(obj.lastDirection,'up')
        		% update by midpoint formula
        		obj.current_rho_bound = (obj.rho_bound_history(end-obj.lag)...
        			+ obj.current_rho_bound)/2;
        		obj.lag = obj.lag + 1;
        	end

        	obj.rho_bound_history(end+1) = obj.current_rho_bound;
        	obj.lastDirection = 'down';
        end

		function rho_update = get_new_rho(obj)
			rho_update = obj.current_rho_bound;
		end
	end

end
