classdef Checks
	% A collection of assertions, mostly for checking
	% numeric arrays.

	methods (Static)
		function has_shape(caller, variable, shape)
			% Throws if variable is not a real array
			% with the required shape.

			if ~isreal(variable)
				error(create_error_struct(caller, "NotReal",...
					"Input is not an array of real numbers"))
			end

            for k = 1:length(shape)
            	if size(variable, k) ~= shape(k)
            		msg = sprintf("Expected length %d along dimension %d%");
            		error(create_error_struct(caller, "InvalidShape", msg))
            	end
            end
		end

		function have_same_shape(caller, varargin)
			% Checks that all inputs are arrays of the same shape.

			if nargin < 3
				return
			end

			for k = 1:nargin
				if ~(isnumeric(varargin{k}) || iscell(varargin{k}))
					error(create_error_struct(caller, "InvalidType",...
						"One or more of the inputs is not an array"))
				end
			end

			shape = size(varargin{1});
			for k = 2:nargin
				if ~isequal(shape, size(varargin{k}))
					error(create_error_struct(caller, "InvalidShape",...
						"Inputs do not have the same shape"))
				end
			end
		end

		function is_square_sparse_matrix(caller, variable, n)
			% Throws if 'variable' is not a real, sparse,
			% nonempty square matrix.

			import HACTLib.Checks

			Checks.is_square_matrix(caller, variable);
			if ~issparse(variable)
				error(create_error_struct(caller, "NotSparse",...
					"Matrix must be sparse"))
			end

			if nargin == 2
				if size(variable, 1) ~= n
					error(create_error_struct(caller, "InvalidShape",...
						"Size of matrix does not match declared size"))
				end
			end
		end

		function is_square_matrix(caller, variable)
			% Throws if 'variable' is not a real, nonempty,
			% square matrix.

			dim1dim2_equal = (size(variable, 1) == size(variable, 2));
			if ~(ismatrix(variable) && dim1dim2_equal)
				error(create_error_struct(caller, "InvalidShape",...
					"Input array is not a square matrix"))
			end

			if size(variable, 1) == 0
				error(create_error_struct(caller, "Empty",...
					"At least one input is empty"))
			elseif ~isreal(variable)
				error(create_error_struct(caller, "NotReal",...
					"Input is complex and must be real"))
			end
		end

		function is_logical(caller, variable)
			if ~islogical(variable)
				msg = "Input must be true/false";
				error(create_error_struct(caller, 'InvalidType', msg));
			end
		end

		function is_real_integer(caller, variable)
			msg = "Input must be real";

			if ~isreal(variable)
				error(create_error_struct(caller, 'NotReal', msg))
			elseif variable ~= round(variable)
				error(create_error_struct(caller, 'InvalidType',...
					"Input is not an integer"))
			end
		end

		function is_integer(caller, variable)
			msg = "Input must be an integer";
			assert(isscalar(variable), msg);
			assert(round(variable == variable), msg);
		end

		function has_attributes(caller, object, props)
			% Checks 'object' to verify that it has at least
			% the required attributes.
		
			assert(iscell(props),...
				"Second arg to has_properties must be a cell array")

			for i = 1:numel(props)
				if ~isprop(object, props{i})
					msg = sprintf("Required attribute '%s' not found", props{i});
					error(msg)
				end
			end
		end
	end
end

function err_struct = create_error_struct(caller, errortype, msg)
	err_struct.identifier = strcat('HACTLib:', caller, ':', errortype);
	err_struct.message = msg;
end