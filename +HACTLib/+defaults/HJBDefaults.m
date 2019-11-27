classdef HJBDefaults
	% Stores the default values for HJBOptions
	properties (Constant)
		implicit = false;
		delta = 1e5;
		HIS_maxiters = 0;
		HIS_tol = 1e-5;
		HIS_start = 2;
	end
end