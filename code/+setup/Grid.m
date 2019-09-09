classdef Grid < handle
    % This class stores the asset and optionally the consumption grids
    
	properties (SetAccess=protected)
		b = struct(); % liquid asset grid
		a = struct(); % illiquid asset grid
		c = struct(); % consumption grid
		z = struct(); % dimension for preference/other heterogeneity

        trapezoidal = struct(); % from trapezoidal rule

        gtype; % HJB or KFE

        ny;
        nz;

        nb;
        nb_pos;
        nb_neg;

        % with pref heterog, rho values are rho + rho_grid
        rho_grid = 0; % e.g. [-0.001,0,0.001]
        rhos; % equals rho + rho_grid

        loc0b;
	end

	methods
	    function obj = Grid(params,ny,gtype)
	    	% pass gtype = 'HJB' or 'KFE'
	    	obj.gtype = gtype;
	    	obj.ny = ny;
            obj.nz = params.nz;

	    	if strcmp(gtype,'KFE')
	    		% Use different-sized KFE grid
                obj.nb = params.nb_KFE;
                obj.nb_neg = params.nb_neg_KFE;
                obj.nb_pos = params.nb_pos_KFE;
	    	elseif strcmp(gtype,'HJB')
                obj.nb = params.nb;
                obj.nb_neg = params.nb_neg;
                obj.nb_pos = params.nb_pos;
            else
	    		error('invalid input')
	    	end
	    end

	    function [grid_vec,grid_wide,grid_mat] = create_generic_grid(...
	    	obj,params,dim,curv,gmin,gmax)

			grid_vec = linspace(0,1,dim)';
			grid_vec = grid_vec.^(1/curv);
			grid_vec = gmin + (gmax - gmin) * grid_vec;
			
			for ia = 1:dim-1
		        if grid_vec(ia+1) - grid_vec(ia) < params.min_grid_spacing
		            grid_vec(ia+1) = grid_vec(ia) + params.min_grid_spacing;
		        else
		            break
		        end
            end

            grid_wide = reshape(grid_vec,[1 dim 1 1]);
			grid_mat = repmat(grid_wide,[obj.nb,1,obj.nz,obj.ny]);
	    end

	    function [bgrid_vec,bgrid_mat] = create_bgrid(obj,params,dim2)
	    	% creates the liquid asset grid
	    	%
	    	% dim2 is the second dimension

	    	% positive part
			bgridpos = linspace(0,1,obj.nb_pos)';
			bgridpos = bgridpos.^(1/params.b_gcurv_pos);
			bgridpos = params.b_soft_constraint + (params.bmax - params.b_soft_constraint) * bgridpos;
			
			for ib = 1:obj.nb_pos
		        if bgridpos(ib+1) - bgridpos(ib) < params.min_grid_spacing
		            bgridpos(ib+1) = bgridpos(ib) + params.min_grid_spacing;
		        else
		            break
		        end
		    end

			% negative part
			if obj.nb_neg > 0
				nb_neg1 = ceil(obj.nb_neg/2) + 1;
				nb_neg2 = obj.nb_neg - nb_neg1 + 2;
				mid_neg = (params.b_soft_constraint + params.bmin) / 2;

				% part of grid close to borrowing limit
				bgridneg1 = linspace(0,1,nb_neg1)';
				bgridneg1 = bgridneg1.^(1/params.b_gcurv_neg);
				bgridneg1 = params.bmin + ...
					(mid_neg - params.bmin) * bgridneg1;
		        
		        % check if min_grid_spacing is too large
		        chunk = (params.b_soft_constraint - params.bmin) / obj.nb_neg;
		        if chunk < params.min_grid_spacing
		            neg_min_spacing = chunk / 2;
		        else
		            neg_min_spacing = params.min_grid_spacing;
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
				bgridneg2 = linspace(1,0,nb_neg2)';
				bgridneg2 = bgridneg2.^(1/params.b_gcurv_neg);
		        bgridneg2 = 1-bgridneg2;
		        bgridneg2 = mid_neg + ...
			    	(params.b_soft_constraint-mid_neg) * bgridneg2;
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
		        
		        bgridneg = [bgridneg1;bgridneg2];

			    bgrid_vec = [bgridneg;bgridpos];
			else
			    bgrid_vec = bgridpos;
			end
			bgrid_mat = repmat(bgrid_vec,[1,dim2,obj.nz,obj.ny]);
		    
		    assert(all(diff(bgrid_vec)>0),'bgrid not strictly increasing')
	    end

	    function [dF,dB,d_tilde] = finite_diff(obj,variable,otherdim)
			% This function finds forward and backward differences, with 
			% variable = 'a','b', or 'c'

			if strcmp(variable,'a')
				input_grid = obj.a.vec;
				repvec = [1 otherdim];
			elseif strcmp(variable,'b')
				input_grid = obj.b.vec;
				repvec = [1 otherdim];
			elseif strcmp(variable,'c')
				input_grid = obj.c.vec;
				repvec = [1 otherdim];
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
		    
		    if strcmp(variable,'a') || strcmp(variable,'c')
		    	% should be in second dimension
		        dF = permute(dF,[2 1]);
		        dB = permute(dB,[2 1]);
		    end
		end

		function add_zgrid(obj, z_vector, dim2)
			obj.z.vec = (1:numel(z_vector))';
			obj.z.wide = reshape(z_vector,[1 1 obj.nz 1]);
			obj.z.matrix = repmat(obj.z.wide,[obj.nb dim2 1 obj.ny]);
		end
    end
end