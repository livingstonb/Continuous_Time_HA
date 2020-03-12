function out = upwind_consumption(net_income_liq_hourly, Vb_fd, direction,...
	prefs, rho_mat, hours_fn)

	endogenous_labor = (nargin > 5);

	if ~ismember(direction, {'B', 'F'})
		error('HACTLib:upwind_consumption:InvalidEntry',...
			'Upwind direction must be B or F')
	end

	nb = size(Vb_fd, 1);

    % Hours worked
    out.hours = hours_fn(Vb_fd);
    out.hours = min(out.hours, 1);

    net_income_liq = net_income_liq_hourly(out.hours);

	out.c = prefs.u1inv(Vb_fd ./ rho_mat);

	if strcmp(direction, 'B')
		out.c(1,:,:,:) = net_income_liq(1,:,:,:);
	else
    	out.c(nb,:,:,:) = 0;
    end

    out.s = net_income_liq - out.c;

    % Impose a state constraint to improve stability
    if strcmp(direction, 'B')
    	out.s(1,:,:,:) = 0; 
    else
    	out.s(nb,:,:,:) = 0;
    end

    out.H = rho_mat .* prefs.u(out.c) + Vb_fd .* out.s ...
        - prefs.hrs_u(out.hours);

    if strcmp(direction, 'F')
        out.H(nb,:,:,:) = -1e12;
    end
end