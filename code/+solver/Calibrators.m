classdef Calibrators
	methods(Static)
        function f = rho_calibrator(runopts, p)
            f = @(x) rho(x, runopts, p);
        end
        
		function f = ra_rho_calibrator(runopts, p)
			f = @(x) ra_rho(x, runopts, p);
		end

		function f = rb_ra_calibrator(runopts, p)
			f = @(x) rb_ra(x, runopts, p);
        end
        
        function x0 = rb_ra_get_initial(p, initial)
            x0 = rb_ra_initial(p, initial);
        end
	end
end

function y = rho(x, runopts, p)
    new_rho = 0.05 * abs(x) / (1 + abs(x));
    p.set("rho", new_rho);
    
    stats = main(runopts, p);
    y = stats.totw - p.targetAY;
    
    fprintf("For rho = %f:\n", p.rho);
    fprintf("Total wealth = %f\n", stats.totw);
end

function y = ra_rho(x, runopts, p)
	new_rho = 0.15 * abs(x(1)) / (1 + abs(x(1)));
    new_ra = p.r_b + 0.1 * abs(x(2)) / (1 + abs(x(2)));
    
    % Set new discount rate
    p.set("rho", new_rho);

	% Set new illiquid returns
	p.set("r_a", new_ra);

    % Solve model
	stats = main(runopts, p);

	% Compute distance from target
    y = zeros(2, 1);
	y(1) = (stats.liqw - 0.5) ^ 2;
    y(2) = (stats.totw - p.targetAY) ^ 2;

    fprintf("For rho=%f:\n", p.rho);
    fprintf("For r_a = %f:\n", p.r_a);
	fprintf("Liquid wealth = %f\n", stats.liqw);
    fprintf("Total wealth = %f\n", stats.totw);
end

function x0 = rb_ra_initial(p, initial)
    if p.riskaver <= 2
		rb_scale = 0.06;
		ra_scale = 0.1;
	elseif p.riskaver <= 10
		rb_scale = 0.13;
		ra_scale = 0.17;
	else
		rb_scale = 0.15;
		ra_scale = 0.25;
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

function [y, x0] = rb_ra(x, runopts, p)
	% This function solves the model for given values of
	% liquid returns and illiquid returns

	if p.riskaver <= 2
		rb_scale = 0.06;
		ra_scale = 0.1;
	elseif p.riskaver <= 10
		rb_scale = 0.13;
		ra_scale = 0.17;
	else
		rb_scale = 0.15;
		ra_scale = 0.25;
	end

	% Set new values for returns
	new_rb = rb_scale*(x(1))/(1+abs(x(1)));
	new_ra = new_rb + ra_scale * abs(x(2)) / (1 + abs(x(2)));
	p.set("r_b", new_rb);
    p.set("r_a", new_ra);
	
	% Solve model
	stats = main(runopts, p);

	% Compute distance from target
	y = (stats.liqw - 0.5) ^ 2;
	y(2) = (stats.totw - p.targetAY) ^ 2;

    fprintf("For r_b = %f", p.r_b);
	fprintf(" and r_a = %f:\n", p.r_a);
	fprintf("Liquid wealth = %f\n", stats.liqw);
	fprintf("Total wealth = %f\n", stats.totw);
end