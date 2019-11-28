function out = upwind_consumption(net_income_liq, Vb_fd, direction,...
	prefs, hours_fn, labor_disutil)

	endogenous_labor = (nargin > 5);

	if ~ismember(direction, {'B', 'F'})
		error('HACTLib:upwind_consumption:InvalidEntry',...
			'Upwind direction must be B or F')
	end

	nb = size(Vb_fd, 1);

	out.c = prefs.u1inv(Vb_fd);

	if strcmp(direction, 'B')
		out.c(1,:,:,:) = net_income_liq(1,:,:,:);
	else
    	out.c(nb,:,:,:) = 0;
    end

    % Hours worked
    if endogenous_labor
        out.hours = hours_fn(Vb_fd);
        out.hours = min(out.hours, 1);
        out.s = net_income_liq(out.hours) - out.c;
    else
        out.s = net_income_liq - out.c;
    end

    % Impose a state constraint to improve stability
    if strcmp(direction, 'B')
    	out.s(1,:,:,:) = 0; 
    else
    	out.s(nb,:,:,:) = 0;
    end

    out.H = prefs.u(out.c) + Vb_fd .* out.s;

    if strcmp(direction, 'F')
        out.H(nb,:,:,:) = -1e12;
    end

    if endogenous_labor
        out.H = out.H - labor_disutil(out.hours);
    end
end