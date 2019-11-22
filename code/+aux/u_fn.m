function u = u_fn(c, gamma)
    % utility function
    %
    % dimz is the dimension in 'c' along which there is heterogeneity
    % in risk aversion. This is a mandatory argument if numel(gamma) > 1
    %
    % the last argument is for stochastic differential utility, optional

    if numel(gamma) == 1
	    if gamma == 1
	        u = log(c);
	    else
	        u = 1/(1-gamma) * (c.^(1-gamma));
	    end
	else
		% need to check which is the z dimension
		ngamma = numel(gamma);
		dims = size(c);

		zdim = find(dims==ngamma,'last');
		gamma1 = (gamma==1);
		
		u = zeros(size(c));
		switch zdim
			case 1
				u(gamma1,:) = log(c(gamma1,:));

				gammasNeq1 = gamma(~gamma1);
				u(~gamma1,:) = 1. / (1-gammasNeq1) * (c(~gamma1,:) .^ (1-gammasNeq1));
			case 2
				u(:,gamma1,:) = log(c(:,gamma1,:));

				gammasNeq1 = reshape(gamma(~gamma1),1,[]);
				u(:,~gamma1,:) = 1. / (1-gammasNeq1) * (c(:,~gamma1,:) .^ (1-gammasNeq1));
			case 3
				u(:,:,gamma1,:) = log(c(:,:,gamma1,:));

				gammasNeq1 = reshape(gamma(~gamma1),1,1,[]);
				u(:,:,~gamma1,:) = 1. / (1-gammasNeq1) * (c(:,:,~gamma1,:) .^ (1-gammasNeq1));
		end
	end
end