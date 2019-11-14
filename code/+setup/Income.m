classdef Income < handle
    % This class stores income grids and related variables

    % Income process is read from files in the directory p.DirIncomeProcess
    % and normalized so that mean quarterly income = 1/4
    
	properties (SetAccess = protected)
		logy;
		y;
		ytrans;
		ydist;
		ny;
	end

	methods
		function obj = Income(runopts,p,dimsHJB,dimsKFE,norisk)
  	 		% this function creates income grids
  	 		%
  	 		% dimsHJB and dimsKFE should be row vectors, e.g.
  	 		% having the values [nb na nz]

  	 		if norisk
  	 			obj.ydist = 1;
  	 			obj.ytrans = 0;

  	 			obj.y.vec = 1/4;
  	 			obj.logy.vec = log(obj.y.vec);
  	 		else
	            incdir = [runopts.direc p.DirIncomeProcess '/'];
	            logy = load([incdir 'ygrid_combined.txt']);
	            y = exp(logy);
	            
	            obj.ydist = load([incdir 'ydist_combined.txt']);
	            obj.ytrans = load([incdir 'ymarkov_combined.txt']);
	            obj.ydist = aux.stat_dist(obj.ytrans');
	            
				% normalize
				obj.y.vec = y ./ (y' * obj.ydist * 4);
				obj.y.vec = obj.y.vec;
				
				obj.logy.vec = log(obj.y.vec);
	        end

	        obj.ny = numel(obj.y.vec);
            obj.y.wide = reshape(obj.y.vec,[1 1 1 obj.ny]);
		    obj.y.matrix = repmat(obj.y.wide,[dimsHJB 1]);
            obj.y.matrixKFE = repmat(obj.y.wide,[dimsKFE 1]);
        end
	end
end