classdef CRRA
	% Provides functions for CRRA utility, with or
	% without IES heterogeneity.

	methods (Static)
		function u = utility(c, invies, zdim)
			% CRRA utility function. Heterogeneity in IES
			% is accomodated by passing as the third argument
			% the dimension of IES heterogeneity in 'c'. If it
			% not passed, an attempt will be made to infer it.

			nIES = numel(invies);

			if nIES == 1
				% No IES heterogeneity
				if invies == 1
					u = log(c);
				else
					u = (1 / (1-invies)) * c .^ (1-invies);
				end
			else
				% IES heterogeneity
				invies1 = (invies == 1);
				inviesNot1 = invies(~invies1);

				if nargin == 2
					zdim = look_for_zdim(c, nIES);
				else
					if size(c, zdim) ~= numel(invies)
						error("HACT:crra:InvalidArgument",...
							strcat("invies dimension of first input does ",...
								"not match the size of the second argument"))
					end
				end

				dims_mask = size(c);
				dims_mask(zdim) = 1;

				invies1_arr = @(iz) logical((invies(iz) == 1) * ones(dims_mask));
				invies1_mask = invies1_arr(1);
				for iz = 2:nIES
					invies1_mask = cat(zdim, invies1_mask, invies1_arr(iz));
				end

				u = zeros(size(c));
				u(invies1_mask) = log(c(invies1_mask));
				
				new_shape = ones(1, ndims(c));
				new_shape(zdim) = nIES;
				rep_shape = size(c);
				rep_shape(zdim) = 1;

				invies_reshaped = reshape(invies, new_shape);
				invies_reshaped = repmat(invies_reshaped, rep_shape);
				inviesNot1 = invies_reshaped(~invies1_mask);
				
				u(~invies1_mask) = 1 ./ (1-inviesNot1) .* (c(~invies1_mask) .^ (1-inviesNot1));
			end
		end

		function muc = marginal_utility(c, invies, zdim)
			nIES = numel(invies);

			if nIES > 1
				if nargin == 2
					zdim = look_for_zdim(c, nIES);
				else
					if size(c, zdim) ~= numel(invies)
						error("HACT:crra:InvalidArgument",...
							strcat("invies dimension of first input does ",...
								"not match the size of the second argument"))
					end
				end

				new_shape = ones(1, ndims(c));
				new_shape(zdim) = nIES;
				invies = reshape(invies, new_shape);
			end

			muc = c .^ (-invies); 
		end

		function c = u1inv(v, invies, zdim)
			nIES = numel(invies);

			if nIES > 1
				if nargin == 2
					zdim = look_for_zdim(c, nIES);
				else
					if size(c, zdim) ~= numel(invies)
						error("HACT:crra:InvalidArgument",...
							strcat("invies dimension of first input does ",...
								"not match the size of the second argument"))
					end
				end

				new_shape = ones(1, ndims(c));
				new_shape(zdim) = nIES;
				invies = reshape(invies, new_shape);
			end

			c = v .^ (-1./invies);
		end
	end
end

function zdim = look_for_zdim(c, nIES)
	if ndims(c) == 4
		zdim = 3;
	else
		zdim = find(size(c)==nIES);
		if numel(zdim) > 1
			error("HACT:crra:InvalidArgument",...
				strcat("Couldn't infer z-dimension, ",...
					"try passing it explicitly"))
		elseif isempty(zdim)
			error("HACT:crra:InvalidArgument",...
				strcat("Could not find a dimension of size ",...
					"numel(invies), check shape of input"))
		end
	end
end