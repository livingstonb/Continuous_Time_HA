classdef ParamsConEffort < setup.Params
    
    properties (SetAccess={?setup.Params})
    	ModelType = 'ConEffort';

    	cmin;
    	cmax;
    	c_gcurv;

    	hdef;

    	chi0;
    	chi1;
    	chi2;
    	hbar;

    	penalty1;
    	penalty2;

    	nc;
    	nc_KFE;

    	delta_HJBguess = 1e5;
	end

    methods
		function obj = ParamsConEffort(runopts,params)
			obj = obj@setup.Params(runopts,params);
            obj.nz = 1;

			if runopts.fast == 1
                obj.nb = 16;
                obj.nb_pos = 14;
                obj.nb_KFE = 15;
                obj.nb_pos_KFE = 14;

                obj.nc = 13;
                obj.nc_KFE = 12;
            end
        end

        function obj = reset_rho(obj,newrho)
            obj.rho = newrho;
        end
        
        function print(obj)
        	fprintf('\n\nSelected parameterization %i:\n',num2str(obj.param_index)) 
            fprintf('%s\n\n',obj.name)

            fprintf('Chosen parameters were...\n\n')

            fprintf('\tb_soft_constraint = %f\n',obj.b_soft_constraint)
            fprintf('\tbmin = %f\n',obj.bmin)
			fprintf('\tbmax = %f\n',obj.bmax)
			fprintf('\tb_gcurv_pos = %f\n',obj.b_gcurv_pos)
			fprintf('\tb_gcurv_neg = %f\n',obj.b_gcurv_neg)
			fprintf('\tnb = %i\n',obj.nb)
			fprintf('\tnb_pos = %i\n',obj.nb_pos)
			fprintf('\tnb_KFE = %i\n',obj.nb_KFE)
			fprintf('\tnb_pos_KFE = %i\n',obj.nb_pos_KFE)
			fprintf('\n')
            
            fprintf('\tcmin = %f\n',obj.cmin)
            fprintf('\tcmax = %f\n',obj.cmax)
            fprintf('\tnc = %i\n',obj.nc)
            fprintf('\tnc_KFE = %i\n',obj.nc_KFE)
            fprintf('\tc_gcurv = %f\n',obj.c_gcurv)
            fprintf('\n')

			fprintf('\tBequests = %i\n',obj.Bequests)
			fprintf('\n')
            
            fprintf('\tchi0 = %f\n',obj.chi0)
            fprintf('\tchi1 = %f\n',obj.chi1)
            fprintf('\tchi2 = %f\n',obj.chi2)
            fprintf('\tpenalty1 = %f\n',obj.penalty1)
            fprintf('\tpenalty2 = %f\n',obj.penalty2)

			fprintf('\tr_b = %f\n',obj.r_b)
			fprintf('\tperfectannuities = %f\n',obj.perfectannuities)
			fprintf('\triskaver = %f\n',obj.riskaver)
			fprintf('\tdeathrate = %f\n',obj.deathrate)
            fprintf('\n\n')
        end
	end
end
