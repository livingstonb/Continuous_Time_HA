classdef Params < HACTLib.model_objects.ParamsBase
    % This class stores the parameters of the model and contains methods to
    % adjust them for frequency and other factors.
    %
    % 
    % See the base class, ParamsBase, for default values.

    methods
        function obj = Params(runopts, varargin)

        	import HACTLib.computation.HJBOptions
        	import HACTLib.computation.KFEOptions

        	obj.set_from_structure(varargin{:});

            % Check for other heterogeneity
            if (numel(obj.rho_grid)>1) && (numel(obj.riskaver)>1)
            	error('Cannot have both rho and riskaver heterogeneity')
            else
                obj.nz = max(numel(obj.rho_grid),numel(obj.riskaver));
            end

             obj.kfe_options = KFEOptions(...
				'maxiters', obj.KFE_maxiters,...
				'tol', obj.KFE_tol,...
				'delta', obj.KFE_delta,...
				'iterative', obj.KFE_iterative);

            obj.hjb_options = HJBOptions(...
				'delta', obj.HJB_delta,...
				'implicit', obj.HJB_implicit,...
				'HIS_maxiters', obj.HIS_maxiters,...
				'HIS_tol', obj.HIS_tol,...
				'HIS_start', obj.HIS_start);
            
            obj.param_index = runopts.param_index;
            obj.ComputeMPCS = runopts.ComputeMPCS;
            obj.SimulateMPCS = runopts.SimulateMPCS;
            obj.ComputeMPCS_news = runopts.ComputeMPCS_news;
            obj.SimulateMPCS_news = runopts.SimulateMPCS_news;
            obj.main_dir = runopts.direc;
            obj.temp_dir = runopts.temp;

            obj.rhos = obj.rho + obj.rho_grid;
            obj.riskaver_fulldim = reshape(obj.riskaver,[1 1 numel(obj.riskaver) 1]);

            obj.SimulateMPCS = runopts.SimulateMPCS;
            obj.DealWithSpecialCase = runopts.DealWithSpecialCase;

            % adjust interest rates
            obj.r_b_borr = obj.r_b + obj.borrwedge;

            % adjust if oneasset is selected
            if obj.OneAsset == 1
                obj.na = 2;
                obj.na_KFE = 2;
                obj.nb = 500;
                obj.nb_KFE = 400;
                obj.b_gcurv_pos = 0.2;
                obj.b_gcurv_neg = 0.2;
                obj.a_gcurv = 0.2;
                obj.bmax = 100;
                obj.chi0 = 1e8;
                obj.r_a = obj.r_b;
            end

            if runopts.fast == 1
                obj.nb = 16;
                obj.nb_pos = 16;
                obj.nb_neg = 0;
                obj.nb_KFE = 15;
                obj.nb_pos_KFE = 15;
                obj.nb_neg_KFE = 0;

                % illiquid grid parameters
                if obj.OneAsset == 0
                    obj.na = 14;
                    obj.na_KFE = 13;
                end
                
                obj.n_mpcsim = 100;
                obj.T_mpcsim = 1e3;
            end

            % Set default grid sizes
            if isempty(obj.nb_pos) == 1
                obj.nb_pos = obj.nb;
            end
            if isempty(obj.nb_pos_KFE) == 1
                obj.nb_pos_KFE = obj.nb_KFE;
            end
            
            obj.nb_neg = obj.nb - obj.nb_pos;
            obj.nb_neg_KFE = obj.nb_KFE - obj.nb_pos_KFE;
        end

        function set_from_structure(obj, varargin)
    		parser = inputParser;

    		import HACTLib.model_objects.ParamsBase

    		defaults = ParamsBase();
    		fields = properties(ParamsBase);

    		for k = 1:numel(fields)
    			field = fields{k};
    			addParameter(parser, field, defaults.(field));
    		end
    		parse(parser, varargin{:});

    		for k = 1:numel(fields)
    			field = fields{k};
    			obj.(field) = parser.Results.(field);
    		end
    	end
        
        function obj = set(obj, field, new_val, quiet)
        	% Sets the value of a parameter.
        	%
        	% Inputs
        	% ------
        	%
        	% field : A string containing the parameter
        	%	name.
        	%
        	% new_val : The desired value of the parameter.
        	%
        	% quiet : An optional argument that, when it
        	%	evaluates to true, suppresses printing
        	%	to the screen.

        	field = char(field);

        	KFE_option_passed = false;
        	if numel(field) > 4
        		if strcmp(field(1:3), 'KFE')
        			KFE_option_passed = true;
        		end
        	end

        	if strcmp(field, 'rho')
        		obj.rhos = new_val + obj.rho_grid;
        	elseif KFE_option_passed
        		assert(isprop(obj.kfe_options, field(5:end)), "Invalid KFE option");
        		obj.kfe_options.set(field(5:end), new_val);
    		elseif ~isprop(obj, field)
    			error("Requested field is not an attribute of Params.");
        	end

            obj.(field) = new_val;

            if ~exist('quiet', 'var')
            	quiet = false;
            end

            if ~quiet
            	disp(strcat(field, sprintf(" has been reset to %.9f", new_val)));
            end
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

            if obj.OneAsset == 0
                fprintf('\tna = %i\n',obj.na)
                fprintf('\tna_KFE = %i\n',obj.na_KFE)
                fprintf('\tamin = %f\n',obj.amin)
                fprintf('\tamax = %f\n',obj.amax)
                fprintf('\ta_gcurv = %f\n',obj.a_gcurv)
                fprintf('\n')

                fprintf('\tchi0 = %f\n',obj.chi0)
                fprintf('\tchi1 = %f\n',obj.chi1)
                fprintf('\tchi2 = %f\n',obj.chi2)
                fprintf('\ta_lb = %f\n',obj.a_lb)
                fprintf('\n')
            end

            fprintf('\tBequests = %i\n',obj.Bequests)
            fprintf('\tResetIncomeUponDeath %i\n',obj.ResetIncomeUponDeath)
            fprintf('\n')

            fprintf('\tr_b = %f\n',obj.r_b)
            fprintf('\tr_b_borr = %f\n',obj.r_b_borr)
            fprintf('\tr_a = %f\n',obj.r_a)
            fprintf('\tperfectannuities = %f\n',obj.perfectannuities)
            fprintf('\triskaver = %f\n',obj.riskaver)
            fprintf('\tdeathrate = %f\n',obj.deathrate)
            fprintf('\n\n')
        end
    end
end
