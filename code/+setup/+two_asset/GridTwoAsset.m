classdef GridTwoAsset < setup.Grid
    
    properties (SetAccess=protected)
        na;
    end
    
	methods
		function obj = GridTwoAsset(params,ny,gtype)
			obj = obj@setup.Grid(params,ny,gtype);

			if strcmp(gtype,'KFE')
	    		% use KFE grid
	    		obj.na = params.na_KFE;
		    elseif strcmp(gtype,'HJB')
		    	% use HJB grid
		    	obj.na = params.na;
	    	else
	    		error('invalid input')
            end

            nz = params.nz;

            % -------------- call methods/create grids----------------
		    % create grids
		    [obj.a.vec, obj.a.matrix] = obj.create_generic_grid(...
	    		params,obj.na,params.a_gcurv,params.amin,params.amax);
		    [obj.b.vec, obj.b.matrix] = obj.create_bgrid(params,obj.na);
		   
            % differenced grids
		    [obj.a.dF, obj.a.dB, da_tilde] = obj.finite_diff('a',obj.nb);
		    [obj.b.dF, obj.b.dB, db_tilde] = obj.finite_diff('b',obj.na);

		    % use trapezoidal rule
		    obj.trapezoidal.vec = kron(da_tilde,db_tilde);
		    obj.trapezoidal.matrix = reshape(repmat(obj.trapezoidal.vec,ny*nz,1),obj.nb,obj.na,ny,nz);
            obj.trapezoidal.matrix_switched = reshape(repmat(obj.trapezoidal.vec,ny*nz,1),obj.nb,obj.na,nz,ny);
			obj.trapezoidal.diagm = spdiags(repmat(obj.trapezoidal.vec,ny*nz,1),0,obj.nb*obj.na*ny*nz,obj.nb*obj.na*ny*nz);
		end
	end

end