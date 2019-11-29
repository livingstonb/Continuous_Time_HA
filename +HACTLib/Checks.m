classdef Checks
	% A collection of assertions, mostly for checking
	% numeric arrays.

	methods (Static)
		function has_shape(caller, variable, shape)
			% Throws if variable is not a real array
			% with the required shape.

			if ~isreal(variable)
				me = gen_exception("Input is not an array of real numbers",...
						caller, "NotReal");
				throwAsCaller(me);
			end

            for k = 1:length(shape)
            	if size(variable, k) ~= shape(k)
            		me = gen_exception("Expected length %d along dimension %d%",...
						caller, "InvalidShape");
					throwAsCaller(me);
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
					me = gen_exception("One or more of the inputs is not an array",...
						caller, "InvalidType");
					throwAsCaller(me);
				end
			end

			shape = size(varargin{1});
			for k = 2:nargin
				if ~isequal(shape, size(varargin{k}))
					me = gen_exception("Inputs do not have the same shape",...
						caller, "InvalidShape");
					throwAsCaller(me);
				end
			end
		end

		function is_square_sparse_matrix(caller, variable, n)
			% Throws if 'variable' is not a real, sparse,
			% nonempty square matrix.

			import HACTLib.Checks

			Checks.is_square_matrix(caller, variable);
			if ~issparse(variable)
				me = gen_exception("Matrix must be sparse",...
						caller, "NotSparse");
				throwAsCaller(me);
			end

			if nargin == 2
				if size(variable, 1) ~= n
					me = gen_exception("Size of matrix does not match declared size",...
						caller, "InvalidShape");
					throwAsCaller(me);
				end
			end
		end

		function is_square_matrix(caller, variable)
			% Throws if 'variable' is not a real, nonempty,
			% square matrix.

			dim1dim2_equal = (size(variable, 1) == size(variable, 2));
			if ~(ismatrix(variable) && dim1dim2_equal)
				me = gen_exception("Input array is not a square matrix", caller,...
					"InvalidShape");
				throwAsCaller(me);
			end

			if size(variable, 1) == 0
				me = gen_exception("At least one input is empty", caller,...
					"Empty");
				throwAsCaller(me);
			elseif ~isreal(variable)
				me = gen_exception("Input must be real", caller,...
					"InvalidType");
				throwAsCaller(me);
			end
		end

		function is_logical(caller, variable)
			if ~islogical(variable)
				me = gen_exception("Input must be true/false", caller,...
					"InvalidType");
				throwAsCaller(me);
			end
		end

		function is_integer(caller, variable)
			if ~isreal(variable)
				me = gen_exception("Input must be real", caller,...
					"InvalidType");
				throwAsCaller(me);
			elseif variable ~= round(variable)
				me = gen_exception("Input is not an integer", caller,...
					"InvalidType");
				throwAsCaller(me);
			end
		end

		function has_attributes(caller, object, props)
			% Checks 'object' to verify that it has at least
			% the required attributes.
		
			if ~iscell(props)
				me = gen_exception("'Attributes' is not a cell array", caller,...
					"InvalidType");
				throwAsCaller(me);
			end

			for i = 1:numel(props)
				if (~isprop(object, props{i}) && ~isfield(object, props{i}))
					msg = sprintf("Required attribute '%s' not found", props{i});
					me = gen_exception(msg, caller, "InvalidType");
					throwAsCaller(me);
				end
			end
		end
	end
end

function err_struct = create_error_struct(caller, errortype, msg)
	err_struct.identifier = strcat('HACTLib:', caller, ':', errortype);
	err_struct.message = msg;
end

function msgId = gen_msgId(caller, varargin)
	msgId = "HACTLib";
	for k = 1:nargin-1
		msgId = cat(msgId, ":", varargin{k});
	end
end

function me = gen_exception(msg, varargin)
	msgId = "HACTLib";
	for k = 1:nargin-1
		msgId = cat(msgId, ":", varargin{k});
	end

	me = MException(msgId, msg);
end