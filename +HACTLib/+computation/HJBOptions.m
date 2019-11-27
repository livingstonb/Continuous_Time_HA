classdef HJBOptions < handle
	% Class used for declaring the options passed to KFESolver.

	properties (SetAccess=private)
		% Set to true for fully implicit updating.
		implicit (1,1) logical...
			= HACTLib.defaults.HJBDefaults.implicit;

		% Step size.
		delta (1,1) double {mustBePositive}...
			= HACTLib.defaults.HJBDefaults.delta;

		% Max number of iterations for the Howard Improvement tep.
		HIS_maxiters (1,1) uint16...
			= HACTLib.defaults.HJBDefaults.HIS_maxiters;

		% Tolerance for the Howard Improvement Step.
		HIS_tol (1,1) double...
			= HACTLib.defaults.HJBDefaults.HIS_tol;

		% Number of HJB iterations before startin the Howard
		% Improvement step.
		HIS_start (1,1) uint16...
			= HACTLib.defaults.HJBDefaults.HIS_start;
	end

	methods
		function obj = HJBOptions(varargin)
			% Creates an HJBOptions object via kew
			options = HACTLib.aux.parse_keyvalue_pairs(...
				HACTLib.defaults.HJBDefaults, varargin{:});
			obj.set_from_struct(options);
		end

		function set_from_struct(obj, options)
			for name = fieldnames(options)'
				obj.(name{1}) = options.(name{1});
			end
		end
	end

end