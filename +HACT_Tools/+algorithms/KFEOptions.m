classdef KFEOptions < handle
	% Class used for declaring the options passed to KFESolver.

	properties (Constant)
		% Default values
		defaults = struct(...
						'iterative', true,...
						'tol', 1e-8,...
						'delta', 1e5,...
						'maxiters', 1e4,...
						'intermediate_check', true...
					);
	end

	properties (SetAccess=private)
		% True if using iterative procedure, false
		% otherwise.
		iterative (1,1) logical = obj.defaults.iterative;

		% Step size for the iterative procedure.
		delta (1,1) double {mustBePositive} = obj.defaults.delta;

		% Convergence tolerance for the iterative
		% procedure.
		tol (1,1) double {mustBePositive} = obj.defaults.tol;

		% Maximum number of iterations for the
		% iterative procedure.
		maxiters (1,1) double {mustBePositive} = obj.defaults.tol;

		% If true, performs a check after a number
		% of iterations. If convergence looks unlikely,
		% an error is thrown.
		intermediate_check (1,1) logical = obj.defaults.intermediate_check;
	end
	methods
		function obj = KFEOptions(varargin)
			% Creates an HJBOptions object via kew
			options = HACT_Tools.aux.parse_keyvalue_pairs(...
				obj.defaults, varargin{:});
			obj.set_from_struct(options);
		end

		function set_from_struct(options)
			for name = fieldnames(options)'
				obj.(name{1}) = options.(name{1});
			end
		end
	end

end