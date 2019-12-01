classdef Grid < handle

	properties (SetAccess=protected)
		vec;
		broadcast;
		array;
	end

	methods
		function obj = Grid(varargin)
		end

		function auto_construct_all(min_val, max_val, n, varargin)
			options = parse_inputs(varargin{:});
				
			obj.auto_construct_vec(min_val, max_val, n, varargin{:});

            if ~isempty(options.dim) && ~isempty(options.state_space)
            	reshape_vec = ones(numel(1, options.state_space));
            	reshape_vec(options.dim) = n;
          		obj.broadcast = reshape(grid_vec, reshape_vec);

          		reshape_vec = options.state_space;
          		reshape_vec(options.dim) = 1;
				obj.a.matrix = repmat(obj.broadcast, reshape_vec);
			end

			assert(all(diff(obj.vec)>0), 'agrid not strictly increasing')
		end

		function auto_construct_vec(min_val, max_val, n, varargin)
			grid_vec = linspace(0, 1, n)';
			grid_vec = grid_vec .^ (1/options.curvature);
			grid_vec = min_val + (max_val - min_val) * grid_vec;
			
			ii = 0;
			while (options.min_spacing > 0) && (ii < n-1)
				ii = ii + 1;
				if grid_vec(ii) - grid_vec(ii) < options.min_spacing
			            grid_vec(ii) = grid_vec(ii) + options.min_spacing;
			        else
			            break
			        end
	            end
            end
           
            obj.vec = grid_vec(:);
        end

		function set
		end


	end

	methods (Access=protected)
		function options = parse_inputs(varargin)
			parser = inputParser;
			addParameters(parser, 'curvature', 1);
			addParameters(parser, 'min_spacing', 0);
			addParameters(parser, 'dim', []);
			addParameters(parser, 'state_space', []);
			parse(parser, varargin{:});
			options = parser.Results;
		end

		function grid_neg = construct_negative_section(min_val, soft_constraint, options)
			n_neg1 = ceil(obj.nb_neg/2) + 1;
			n_neg2 = obj.nb_neg - nb_neg1 + 2;
			mid_neg = (obj.p.b_soft_constraint + obj.p.bmin) / 2;

			% Section of grid close to the borrowing limit
			grid_neg1 = linspace(0, 1, n_neg1)';
			grid_neg1 = grid_neg1.^(1/options.curvature);
			grid_neg1 = obj.p.bmin + ...
				(mid_neg - obj.p.bmin) * grid_neg1;
	        
	        % check if min_grid_spacing is too large
	        chunk = (obj.p.b_soft_constraint - obj.p.bmin) / obj.nb_neg;
	        if chunk < obj.p.min_grid_spacing
	            neg_min_spacing = chunk / 2;
	        else
	            neg_min_spacing = obj.p.min_grid_spacing;
	        end
	        
	        % enforce very small minimum grid spacing
	        for ib = 1:nb_neg1-1
	            if bgridneg1(ib+1) - bgridneg1(ib) < neg_min_spacing
	                bgridneg1(ib+1) = bgridneg1(ib) + neg_min_spacing;
	            else
	                break
	            end
	        end

			% part of grid close to soft constraint
			bgridneg2 = linspace(1, 0, nb_neg2)';
			bgridneg2 = bgridneg2 .^ (1/obj.p.b_gcurv_neg);
	        bgridneg2 = 1-bgridneg2;
	        bgridneg2 = mid_neg + ...
		    	(obj.p.b_soft_constraint - mid_neg) * bgridneg2;
		   	bgridneg2(1) = []; % remove midpoint
		    bgridneg2(end) = []; % remove 0

	        % enforce minimum grid spacing
		    if 0 - bgridneg2(nb_neg2-2) < neg_min_spacing
		        bgridneg2(nb_neg2-2) = - neg_min_spacing;
		    end

		    for ib = nb_neg2-2:-1:2
		        tooClose = bgridneg2(ib) - bgridneg2(ib-1) < neg_min_spacing;
		        invalid = bgridneg2(ib) < bgridneg2(ib-1);
		        if tooClose || invalid
		            bgridneg2(ib-1) = bgridneg2(ib) - neg_min_spacing;
		        else
		            break
	            end
	        end
	        
	        bgridneg = [bgridneg1; bgridneg2];
		    end
		end
	end

end