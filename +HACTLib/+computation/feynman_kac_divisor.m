function B = feynman_kac_divisor(p, income, stepsize, A, invert)
	% Computes the (k,k')-income blocks down the main diagonal for
	% left-hand-size of the Feynman-Kac equation.
	%
	% Required Inputs
	% ---------------
	% p : An object with the following attributes:
	%	nb_KFE, na_KFE, nz
	%	- Grid sizes.
	%
	%	deathrate
	%	- The Poisson death rate.
	%
	% income : An object with the following attributes:
	%	ytrans
	%	- The square income transition matrix. Should
	%	  have row sums of zero.
	%
	%	ny
	%	- The number of income states.
	%
	%	stepsize
	%	- The time step for iteration. Smaller values increase
	%	  the accuracy of the results but will increase
	%	  computational time.
	%
	% Optional Inputs
	% ---------------
	% invert : Boolean indicator for the return type. If
	%	true, then an object B is returned s.t. 
	%	B{k} * A is roughly equivalent to the backslash
	%	operator. If false, B is the array on the
	%	left-hand-side of the Feynman-Kac equation.

	if nargin == 4
		invert = false;
	end

	B = cell(income.ny, 1);
	states_per_income = p.nb_KFE * p.na_KFE * p.nz;
	for k = 1:income.ny
		ind1 = 1 + states_per_income * (k-1);
    	ind2 = states_per_income * k;

        B{k} = speye(states_per_income) * (...
        			1 / stepsize + p.deathrate - income.ytrans(k,k)...
                	) - A(ind1:ind2,ind1:ind2);

        if invert
        	B{k} = inverse(B{k});
        end
	end
end