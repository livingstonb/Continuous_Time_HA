classdef GridConEffort < setup.Grid
    
    properties (SetAccess=protected)
        nc;
    end
    
	methods
		function obj = GridConEffort(params,ny,gtype)
			obj = obj@setup.Grid(params,ny,gtype);

			if strcmp(gtype,'KFE')
	    		% use KFE grid
	    		obj.nc= params.nc_KFE;
		    elseif strcmp(gtype,'HJB')
		    	% use HJB grid
		    	obj.nc = params.nc;
	    	else
	    		error('invalid input')
            end

            nz = params.nz;

            % -------------- call methods/create grids----------------
		    % create grids
		    [obj.c.vec, obj.c.wide, obj.c.matrix] = obj.create_generic_grid(...
	    		params,obj.nc,params.c_gcurv,params.cmin,params.cmax);
		    [obj.b.vec, obj.b.matrix] = obj.create_bgrid(params,obj.nc);
		   
            % differenced grids
		    [obj.c.dF, obj.c.dB, dc_tilde] = obj.finite_diff('c',obj.nb);
		    [obj.b.dF, obj.b.dB, db_tilde] = obj.finite_diff('b',obj.nc);

		    % use trapezoidal rule
		    obj.trapezoidal.vec = kron(dc_tilde,db_tilde);
		    obj.trapezoidal.matrix = reshape(repmat(obj.trapezoidal.vec,ny*nz,1),obj.nb,obj.nc,nz,ny);
			obj.trapezoidal.diagm = spdiags(repmat(obj.trapezoidal.vec,ny*nz,1),0,obj.nb*obj.nc*ny*nz,obj.nb*obj.nc*ny*nz);

			% find location of zero assets in b vector (for death without bequests)
			obj.loc0b = find(obj.b.vec==0);
		end
	end

end