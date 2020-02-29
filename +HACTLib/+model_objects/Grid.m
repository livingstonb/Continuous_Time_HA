classdef Grid < handle
    % This class creates and stores the asset grids.
    % There are two recommended ways to use this class:
    %
    %	(1) Call the auto_construct() after instantiating
    %		a Grid object.
    %
    %	(2) Construct the b and a grids manually as vectors
    %		and pass them to create_agrids() and create_bgrids(),
    %		then call the generate_variables() method.
   
	properties (SetAccess=protected)
		%	Model parameters. Must contain at least the
		%	following fields:
    	%
    	%		nb (nb_KFE if gtype = 'KFE')
    	%		-  Number of grid points for the liquid asset.
    	%
    	%		nb (nb_KFE if gtype = 'KFE')
    	%		- Number of grid points for the illiquid asset.
    	%		nz 
    	%		- Number of grid points in extra heterogeneity
    	%		  dimension, typically 1.
    	%
    	%	and if auto_construct() will be used, the following
    	%	fields are also required:
    	%	
    	%		nb_pos / nb_pos_KFE
    	%		nb_neg / nb_neg_KFE
        p;
        
        % The liquid asset grids.
		b = struct();

		% The illiquid asset grids.
		a = struct();

		% Total wealth.
		w = struct();

		% Dimension for preference/other heterogeneity.
		z = struct();

		% Number of points on the illiquid asset grid.
		na;

		% Number of points on the liquid asset grid.
        nb;

        % Number of points on the non-negative section of
        % the liquid asset grid.
        nb_pos;

        % Number of points on the negative section of
        % the liquid asset grid.
        nb_neg;

        % Linear index for b = 0 along the second dimension.
        loc0b;

		% Linear index for a = 0 along the second dimension.
        loc0a;

        % Linear index for b = 0 and a = 0 in a stacked b/a
        % vector.
        % i.e. repmat(b, na) == 0 and kron(agrid, ones(nb, 1)) == 0
        loc0b0a;
	    

        % Grid deltas for trapezoidal integration.
        da_tilde;
        db_tilde;
        trapezoidal = struct();

        % String indicating which grid is being created, 'HJB'
        % or 'KFE'.
        gtype;

        % Number of points on the income grid.
        ny;

        % Number of states in the extra dimension of
		% heterogeneity.
        nz;

        % Grid for the values of the discount factor, rho,
        % before adding the rho.
        rho_grid = 0; % e.g. [-0.001,0,0.001]

        % A grid of rho values. In particular, the sum of
        % rho and rho_grid, a vector.
        rhos;

        % Indicator of whether or not the all of the grids
        % have been initialized.
        fully_initialized = false;
	end

	methods
	    function obj = Grid(p, ny, gtype)
	    	% Class constructor.
	    	%
	    	% Required Parameters
	    	% -------------------
	    	% p : An object containing the model parameters. See
	    	%	the help documentation of 'p' for more information.
	    	%
	    	% ny : The number of points on the income grid.
	    	%
	    	% Optional Parameters
	    	% -------------------
	    	% gtype : Type of grid, specify 'KFE' to construct the grids for the KFE.
	    	%	The user may want different grid sizes for the HJB and KFE.

            obj.p = p;

	    	% pass gtype = 'HJB' or 'KFE'
	    	obj.gtype = gtype;
	    	obj.ny = ny;
            obj.nz = p.nz;

	    	if strcmp(gtype,'KFE')
	    		% Use different-sized KFE grid
                obj.nb = p.nb_KFE;

                if isprop(p, 'nb_neg_KFE')
                	obj.nb_neg = p.nb_neg_KFE;
                end

                if isprop(p, 'nb_pos_KFE')
                	obj.nb_pos = p.nb_pos_KFE;
                end
                
                obj.na = p.na_KFE;
	    	else
                obj.nb = p.nb;

                if isprop(p, 'nb_neg')
                	obj.nb_neg = p.nb_neg;
                end

                if isprop(p, 'nb_pos')
                	obj.nb_pos = p.nb_pos;
                end

                obj.na = p.na;
            end
        end
        
        function obj = auto_construct(obj)
            % Automatically constructs all the
            % necessary grids from the parameters
            % used to initialize the Grid object.
            
            obj.create_agrids();
            obj.create_bgrids();
            obj.generate_variables();
            obj.fully_initialized = true;
        end

        function obj = generate_variables(obj)
        	assert(~isempty(obj.a), "agrid has not been constructed")
        	assert(~isempty(obj.b), "bgrid has not been constructed")

        	obj.finite_diff('a');
		    obj.finite_diff('b');

			obj.loc0b = find(obj.b.vec==0);
			obj.loc0a = find(obj.a.vec==0);

			% location of b = 0 and a = 0 in (b,a)-space
			bLong = repmat(obj.b.vec, obj.na, 1);
			aLong = kron(obj.a.vec, ones(obj.nb, 1));
			obj.loc0b0a = find((bLong==0) & (aLong==0));
            
            obj.construct_trapezoidal_grid();
            obj.create_zgrids();
            obj.fully_initialized = true;
        end

	    function create_agrids(obj, grid_vec)
	    	% Creates the illiquid asset grids.
	    	%
	    	% Optional Parameters
	    	% ----------
	    	% grid_vec : Optionally, the user can supply the
	    	%	illiquid asset grid as a vector, and the
	    	%	other variables will be automatically created.
	    	%	This vector must include 0.
	    	%
	    	% Modifies
	    	% -------
	    	% obj.a.vec :  A column vector grid.
	    	%
	    	% obj.a.wide : An array grid of shape (1, na, 1, 1).
	    	%
	    	% obj.a.matrix : An array the full size of the model, (nb, na, nz, ny)

	    	if nargin == 1
				grid_vec = linspace(0, 1, obj.na)';
				grid_vec = grid_vec .^ (1/obj.p.a_gcurv);
				grid_vec = obj.p.amin + (obj.p.amax - obj.p.amin) * grid_vec;
				
				for ia = 1:obj.na-1
			        if grid_vec(ia+1) - grid_vec(ia) < obj.p.min_grid_spacing
			            grid_vec(ia+1) = grid_vec(ia) + obj.p.min_grid_spacing;
			        else
			            break
			        end
	            end
            elseif nargin > 1
	        	assert(isvector(grid_vec), "Input agrid must be a vector")
	        	assert(sum(grid_vec==0) == 1, "Input grid must include 0")
	        end
            obj.a.vec = grid_vec(:);
            obj.a.wide = reshape(grid_vec, [1 obj.na 1 1]);
			obj.a.matrix = repmat(obj.a.wide, [obj.nb, 1, obj.nz, obj.ny]);
			assert(all(diff(obj.a.vec)>0), 'agrid not strictly increasing')
	    end

	    function create_bgrids(obj, bgrid_vec)
	    	% Creates the liquid asset grids.
	    	%
	    	% Optional Parameters
	    	% -------------------
	    	% bgrid_vec : Optionally, the user can supply the
	    	%	liquid asset grid as a vector, and the
	    	%	other variables will be automatically created.
	    	%	This vector must include 0.
	    	%
	    	% Modifies
	    	% --------
	    	% obj.b.vec : A column vector for the liquid asset grid.
	    	%
	    	% obj.b.matrix : An array of the liquid asset grid,
	    	%	shape (nb, na, nz, ny).

	    	msg = ['Number of non-negative and negative liquid grid points',...
	    			'do not sum to total number of liquid grid points'];
	    	assert(obj.nb == obj.nb_pos + obj.nb_neg, msg)
	    
	    	if nargin == 1
		    	% positive part
				bgridpos = linspace(0,1,obj.nb_pos)';
				bgridpos = bgridpos.^(1/obj.p.b_gcurv_pos);
				bgridpos = obj.p.b_soft_constraint + (obj.p.bmax - obj.p.b_soft_constraint) * bgridpos;
				
				for ib = 1:obj.nb_pos
			        if bgridpos(ib+1) - bgridpos(ib) < obj.p.min_grid_spacing
			            bgridpos(ib+1) = bgridpos(ib) + obj.p.min_grid_spacing;
			        else
			            break
			        end
			    end

				% negative part
				if obj.nb_neg > 0
					nb_neg1 = ceil(obj.nb_neg/2) + 1;
					nb_neg2 = obj.nb_neg - nb_neg1 + 2;
					mid_neg = (obj.p.b_soft_constraint + obj.p.bmin) / 2;

					% part of grid close to borrowing limit
					bgridneg1 = linspace(0, 1, nb_neg1)';
					bgridneg1 = bgridneg1.^(1/obj.p.b_gcurv_neg);
					bgridneg1 = obj.p.bmin + ...
						(mid_neg - obj.p.bmin) * bgridneg1;
			        
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

				    obj.b.vec = [bgridneg; bgridpos];
				else
				    obj.b.vec = bgridpos;
				end
			else
				assert(isvector(bgrid_vec), "Input bgrid must be a vector")
				obj.b.vec = bgrid_vec(:);
			end
			obj.b.matrix = repmat(obj.b.vec, [1 obj.na obj.nz obj.ny]);
		    assert(all(diff(obj.b.vec)>0), 'bgrid not strictly increasing')
	    end

	    % function create_wgrids(obj)
	    % 	obj.w.matrix = obj.b.matrix = obj.a.matrix;

	    % end

	    function obj = set(obj, varname, value)
	    	obj.(varname) = value;
	    end
	end

	methods (Access=private)
	    function obj = finite_diff(obj, variable)
			% This function finds forward and backward differences.
			%
			% Parameters
			% ----------
			% variable : either the liquid asset, 'b', or illiquid, 'a'
			%
			% Modifies
			% -------
			% obj.<variable>.dF : the forward first difference of the selected grid
			%
			% obj.<variable>.dB : the backward first difference of the selected grid
			%
			% obj.d<variable>_tilde : the mean of dF and dB, adjusted at the boundaries

			if strcmp(variable,'a')
				input_grid = obj.a.vec;
				repvec = [1 obj.nb];
			elseif strcmp(variable,'b')
				input_grid = obj.b.vec;
				repvec = [1 obj.na];
			else
				error('invalid variable choice')
			end

			% Forward diff
			dF = NaN(size(input_grid));
			dF(1:end-1) = input_grid(2:end) - input_grid(1:end-1);
			dF(end) = dF(end-1);

			% Backward diff
			dB = NaN(size(input_grid));
			dB(2:end) = input_grid(2:end) - input_grid(1:end-1);
			dB(1) = dB(2);

			% Trapezoidal rule
		    d_tilde = (dF + dB)/2;
			d_tilde(1) = dF(1)/2;
			d_tilde(end) = dB(end)/2;
		    
		    dF = repmat(dF,repvec);
		    dB = repmat(dB,repvec);
		    
		    if strcmp(variable,'a')
		    	obj.da_tilde = d_tilde;
		        obj.a.dF = permute(dF,[2 1]);
		        obj.a.dB = permute(dB,[2 1]);
		    else
		    	obj.db_tilde = d_tilde;
		    	obj.b.dF = dF;
		    	obj.b.dB = dB;
		    end
		end

		function construct_trapezoidal_grid(obj)
			% Constructs the grid of db * da.

			obj.trapezoidal.vec = kron(obj.da_tilde,obj.db_tilde);
		    obj.trapezoidal.matrix = reshape(repmat(obj.trapezoidal.vec,obj.ny*obj.nz,1),...
		    	obj.nb,obj.na,obj.nz,obj.ny);
			obj.trapezoidal.diagm = spdiags(repmat(obj.trapezoidal.vec,obj.ny*obj.nz,1), 0,...
				obj.nb*obj.na*obj.ny*obj.nz,obj.nb*obj.na*obj.ny*obj.nz);
		end


		function create_zgrids(obj)
			% Creates an integer index grid for the z-dimension.

			obj.z.vec = (1:obj.nz)';
			obj.z.wide = reshape(obj.z.vec,[1 1 obj.nz 1]);
			obj.z.matrix = repmat(obj.z.wide,[obj.nb obj.na 1 obj.ny]);
        end
	end
end