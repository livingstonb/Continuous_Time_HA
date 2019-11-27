classdef KFEDefaults
	% Stores the default values for KFEOptions
	properties (Constant)
		iterative = true;
		tol = 1e-8;
		delta = 1e5;
		maxiters = 1e4;
		intermediate_check = true;
	end
end