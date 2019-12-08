classdef IncomeSimulator

	properties (SetAccess=protected)
		% Vector of current income indices
		yinds;

		% Sample size
		nsim

		% Income transition matrix
		cum_ytrans

		% Poisson rate of death
		deathrate;
	end

	methods
		function obj = IncomeSimulator(nsim, ytrans, deathrate)
			obj.nsim = nsim;
			obj.cum_ytrans = cumsum(ytrans, 2);

			if nargin == 3
				obj.deathrate = deathrate;
			else
				obj.deathrate = 0;
			end
		end

		function next(obj, incdraws)
			[~, obj.yinds] = max(incdraws...
				<= obj.cum_ytrans(obj.yinds,:), [], 2);
		end

		function simulate_one_period(obj, draws)
			chunksize = 5e2;
			finished = false;
			i1 = 1;
			i2 = min(chunksize,obj.p.n_mpcsim);
			while ~finished
		        [~,obj.yinds(i1:i2)] = max(draws(i1:i2,1)...
		        	<=obj.cum_ytrans(obj.yinds(i1:i2),:), [], 2);

		        if obj.p.Bequests == 0
		        	lived = draws(i1:i2,2) > obj.deathrateSubperiod;
		        	obj.bsim(i1:i2) = lived .* obj.bsim(i1:i2);

		        	if strcmp(obj.dim2Identity,'a')
		        		obj.asim(i1:i2) = lived .* obj.asim(i1:i2);
		        	end
		        end
		        
		        i1 = i2 + 1;
		        i2 = min(i1+chunksize, obj.p.n_mpcsim);
		        
		        if i1 > obj.p.n_mpcsim
		            finished = true;
		        end
			end

	        obj.ysim = obj.income.y.vec(obj.yinds);
            obj.yrep = repmat(obj.ysim, obj.nshocks+1,1);
	    end
	end

end