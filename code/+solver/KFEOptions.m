classdef KFEOptions
	% Class used for declaring the options passed to KFESolver.

	properties (Constant)
		% Property defaults

		default_iterative (1,1) logical = true;

		default_tol (1,1) double {mustBePositive} = 1e-8;

		default_delta (1,1) double {mustBePositive} = 1e5;

		default_maxiters (1,1) double {mustBePositive} = 1e4;

		default_intermediate_check (1,1) logical = true;
	end

	properties (SetAccess=private)
		% True if using iterative procedure, false
		% otherwise.
		iterative (1,1) logical = solver.KFEOptions.default_iterative;

		% Step size for the iterative procedure.
		delta (1,1) double {mustBePositive} = solver.KFEOptions.default_delta;

		% Convergence tolerance for the iterative
		% procedure.
		tol (1,1) double {mustBePositive} = solver.KFEOptions.default_tol;

		% Maximum number of iterations for the
		% iterative procedure.
		maxiters (1,1) double {mustBePositive} = solver.KFEOptions.default_maxiters;

		% If true, performs a check after a number
		% of iterations. If convergence looks unlikely,
		% an error is thrown.
		intermediate_check (1,1) logical = solver.KFEOptions.default_intermediate_check;
	end
	methods
		function obj = KFEOptions(varargin)

			parser = inputParser;
			addParameter(parser, 'delta', obj.default_delta);
			addParameter(parser, 'tol', obj.default_tol);
			addParameter(parser, 'maxiters', obj.default_maxiters);
			addParameter(parser, 'iterative', obj.default_iterative);
			addParameter(parser, 'intermediate_check', obj.default_intermediate_check);

			parse(parser, varargin{:});

			for name = fieldnames(parser.Results)'
				obj.(name{:}) = parser.Results.(name{:});
			end
		end
	end

end