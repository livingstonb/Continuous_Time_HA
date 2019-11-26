function outparams = chi1_chi2_tests(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,...
    
    % chi1s = 0.2:0.2:1;
    chi1s = [0.01,0.05,0.1,0.15];
    chi2s = [0.25,0.5];

    ras_chi2_025_lowchi1 = [0.0218531297928173;0.0202034161839398;0.0195184891122385;0.0190643215902504];

    ras_chi2_05_lowchi1 = [0.0203;0.0183763781021196;0.0169338416193918;0.0158630015792344];
    
    ras_chi2_025 = [0.0187053516216175;0.0177363614686727;0.0170710857984337;0.0165970092805796;0.0162375425622935];
                
    ras_chi2_05 = [0.0150355107927622;0.0130567686762534;0.0119786881432546;0.0112706123306245;0.0107578146068646];
    
    ras_chi2_1 = [0.0103095789849221;0.00820652079041105;0.00758552547233842;0.00835749610168757;0.00973569589474840];
    
    ras_chi2_2 = [0.0085,0.0065,0.0064,0.0063,0.0062];
            
    ii = 1;
    for eps = [0]
    for chi2 = chi2s
        
        chi1Index = 1;
        for chi1 = chi1s
            chi0 = 0;

            params(ii).name = sprintf('contB, chi1_%0.2f, chi2_%0.2f',chi1,chi2); 
            params(ii).OneAsset = 0;
            params(ii).DirIncomeProcess = 'IncomeGrids/continuous_b';
            params(ii).chi2 = chi2;
            params(ii).chi1 = chi1;
            params(ii).chi0 = chi0;
            params(ii).a_lb = 0.25;
            params(ii).rhoL = 0.005;

            if chi2 == 0.25
                params(ii).r_a = ras_chi2_025_lowchi1(chi1Index) + eps;
            elseif chi2 == 0.5
                params(ii).r_a = ras_chi2_05_lowchi1(chi1Index) + eps;
            end
            % j = mod(ii,5);
            % if j == 0
            %     j = 5;
            % end

            % switch chi2
            %     case 0.25
            %         params(ii).r_a = ras_chi2_025(j);
            %     case 0.5
            %         params(ii).r_a = ras_chi2_05(j);
            %     case 1
            %         params(ii).r_a = ras_chi2_1(j);
            %     case 2
            %         params(ii).r_a = ras_chi2_2(j);
            % end

            ii = ii + 1;
            chi1Index = chi1Index + 1;
        end
    end
    end

    %% DO NOT CHANGE BELOW

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACT_Tools.components.Params(runopts,chosen_param);

end
