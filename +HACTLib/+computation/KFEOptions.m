classdef KFEOptions < handle
	% Class used for declaring the options passed to KFESolver.

	properties (SetAccess=private)
		% True if using iterative procedure, false
		% otherwise.
		iterative (1,1) logical...
			= HACTLib.defaults.KFEDefaults.iterative;

		% Step size for the iterative procedure.
		delta (1,1) double...
			= HACTLib.defaults.KFEDefaults.delta;

		% Convergence tolerance for the iterative
		% procedure.
		tol (1,1) double {mustBePositive}...
			= HACTLib.defaults.KFEDefaults.tol;

		% Maximum number of iterations for the
		% iterative procedure.
		maxiters (1,1) double {mustBePositive}...
			= HACTLib.defaults.KFEDefaults.maxiters;

		% If true, performs a check after a number
		% of iterations. If convergence looks unlikely,
		% an error is thrown.
		intermediate_check (1,1) logical...
			= HACTLib.defaults.KFEDefaults.intermediate_check;
	end
	methods
		function obj = KFEOptions(varargin)
			% Creates an HJBOptions object via key-value pairs
			options = HACTLib.aux.parse_keyvalue_pairs(...
				HACTLib.defaults.KFEDefaults, varargin{:});
			obj.set_from_struct(options);
		end

		function set_from_struct(obj, options)
			for name = fieldnames(options)'
				obj.(name{1}) = options.(name{1});
			end
		end
	end

end