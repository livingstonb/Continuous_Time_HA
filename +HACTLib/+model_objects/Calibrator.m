classdef Calibrator
	% Class used in coordination with fsolve to
	% calibrate the model to one or more targets

	properties (SetAccess = protected)
	% objective function for fsolve
	objective;

	% function to get initial conditions for fsolve
	create_initial_condition;
	end
	
	methods
		function obj = Calibrator(runopts, p, method)
			% Parameters
			% ----------
			% runopts : structure containing run options
			%
			% p : a Params object
			%
			% method : a string variable referring to the type of
			%	calibrator requested
			%
			% Returns
			% -------
			% obj : a Calibrator object
			
			if strcmp(method, "rho")
				obj.objective = @(x) rho_calibrator(x, runopts, p);
				obj.create_initial_condition = @(initial_cond) rho_initial(p, initial_cond);
			elseif strcmp(method, "r_a, rho");
				obj.objective = @(x) ra_rho_calibrator(x, runopts, p);
				obj.create_initial_condition = @(initial_cond) ra_rho_initial(p, initial_cond);
			elseif  strcmp(method, "r_b, r_a");
				obj.objective = @(x) rb_ra_calibrator(x, runopts, p);
				obj.create_initial_condition = @(initial_cond) rb_ra_initial(p, initial_cond);
			else
				error("Requested calibrator does not exist");
			end
		end
	end
end

function x0 = rho_initial(p, initial)
	% a function used to create initial conditions for
	% fsolve, out of a desired initial rho

	% Parameters
	% ----------
	% p : a Params object
	%
	% initial : a vector containing the desired initial
	%	value for rho
	%
	% Returns
	% ------
	% x0 : a scalar containing the required input value
	%	for the calibrator to set rho to its
	%	desired value, passed in 'initial'

	x0 = (initial / 0.05) / (1 - initial / 0.05);
end


function y = rho_calibrator(x, runopts, p)
	% objective function for calibration of rho to
	% a mean wealth target

	% Parameters
	% ----------
	% x : a variable to be passed by fsolve; this
	%	variable determines the value of rho
	%	according to a bounded transformation below
	%
	% runopts : a structure of run options
	%
	% p : a Params object
	%
	% Returns
	% -------
	% y : the deviation of actual mean wealth from
	%	target wealth

    new_rho = 0.05 * abs(x) / (1 + abs(x));
    p.set("rho", new_rho);
    
    stats = main(runopts, p);
    y = stats.totw - p.targetAY;
    
    fprintf("For rho = %f:\n", p.rho);
    fprintf("Total wealth = %f\n", stats.totw);
end

function x0 = ra_rho_initial(p, initial)
	% a function used to create initial conditions for
	% fsolve, out of a desired initial rho and r_a

	% Parameters
	% ----------
	% p : a Params object
	%
	% initial : a vector containing the desired initial
	%	values for rho and r_a
	%
	% Returns
	% ------
	% x0 : a vector containing the required input values
	%	for the calibrator to set rho and r_a to their
	%	desired values, passed in 'initial'

    x0 = zeros(2, 1);
   
    x0(1) = ((initial(1)-p.r_b) / 0.25)...
    	/ (1 - (initial(1)-p.r_b) / 0.25);

     if initial(1) > 0
        x0(2) = (initial(2) / 0.35) / (1 - initial(2)/0.35);
    else
        x0(2) = (initial(2) / 0.35) / (1 + initial(2)/0.35);
    end
end


function y = ra_rho_calibrator(x, runopts, p)
	% the objective function when calibrating r_a and rho
	% to match mean wealth to targetAY and mean liquid wealth
	% to 0.5

	% Parameters
	% ----------
	% x : the vector passed from fsolve, to be transformed
	%	to rho and r_a
	%
	% runopts : a structure containing run options
	%
	% p : a Params object
	%
	% Returns
	% -------
	% y : a vector containing the squared deviations of
	%	mean liquid wealth and mean wealth from their
	%	targets

    new_ra = p.r_b + 0.25 * abs(x(1)) / (1 + abs(x(1)));
    new_rho = 0.35 * abs(x(2)) / (1 + abs(x(2)));
    
    % Set new discount rate
    p.set("rho", new_rho);

	% Set new illiquid returns
	p.set("r_a", new_ra);

    % Solve model
	stats = main(runopts, p);

	% Compute distance from target
    y = zeros(2, 1);
	y(1) = stats.liqw - 0.5;
    y(2) = stats.totw - p.targetAY;

    fprintf("For rho=%f:\n", p.rho);
    fprintf("For r_a = %f:\n", p.r_a);
	fprintf("Liquid wealth = %f\n", stats.liqw);
    fprintf("Total wealth = %f\n", stats.totw);
end

function x0 = rb_ra_initial(p, initial)
	% a function used to create initial conditions for
	% fsolve, out of a desired initial r_b and r_a

	% Parameters
	% ----------
	% p : a Params object
	%
	% initial : a vector containing the desired initial
	%	values for r_b and r_a
	%
	% Returns
	% ------
	% x0 : a vector containing the required input values
	%	for the calibrator to set r_b and r_a to their
	%	desired values, passed in 'initial'

    if p.riskaver <= 2
		rb_scale = 0.15;
		ra_scale = 0.15;
	elseif p.riskaver <= 10
		rb_scale = 0.4;
		ra_scale = 0.4;
	else
		rb_scale = 0.5;
		ra_scale = 0.5;
	end
    
     % return initial conditions
    x0 = zeros(2, 1);
    if initial(1) > 0
        x0(1) = (initial(1) / rb_scale) / (1 - initial(1)/rb_scale);
    else
        x0(1) = (initial(1) / rb_scale) / (1 + initial(1)/rb_scale);
    end
    
    x0(2) = ((initial(2)-initial(1)) / ra_scale) / (1 - (initial(2)-initial(1)) / ra_scale);
end

function [y, x0] = rb_ra_calibrator(x, runopts, p)
	% the objective function when calibrating r_a and rho
	% to match mean wealth to targetAY and mean liquid wealth
	% to 0.5

	% Parameters
	% ----------
	% x : the vector passed from fsolve, to be transformed
	%	to r_b and r_a
	%
	% runopts : a structure containing run options
	%
	% p : a Params object
	%
	% Returns
	% -------
	% y : a vector containing the squared deviations of
	%	mean liquid wealth and mean wealth from their
	%	targets

	if p.riskaver <= 2
		rb_scale = 0.15;
		ra_scale = 0.15;
	elseif p.riskaver <= 10
		rb_scale = 0.4;
		ra_scale = 0.4;
	else
		rb_scale = 0.5;
		ra_scale = 0.5;
	end

	% Set new values for returns
	new_rb = rb_scale*(x(1))/(1+abs(x(1)));
	new_ra = new_rb + ra_scale * abs(x(2)) / (1 + abs(x(2)));
	p.set("r_b", new_rb);
    p.set("r_a", new_ra);
	
	% Solve model
	stats = main(runopts, p);

	% Compute distance from target
	y = stats.liqw - 0.5;
	y(2) = stats.totw - p.targetAY;

    fprintf("For r_b = %f", p.r_b);
	fprintf(" and r_a = %f:\n", p.r_a);
	fprintf("Liquid wealth = %f\n", stats.liqw);
	fprintf("Total wealth = %f\n", stats.totw);
end