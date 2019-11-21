function outparams = SDU_tests(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,.

    % rho for calibration based on baseline: 0.033940
    % rho for calibration based on riskaver = 5, sigma_r = 0.1: 0.114248

%     rho_ies1 = 0.0339400;
%     rho_ies1_5 = 0.039426;

    rho_ies1 = 0.1;
%     rho_ies1 = 0.274541;
    rho_ies1_5 = 0.223813;

    %%--------------------------------------------------------------
    % BASELINE
    % --------------------------------------------------------------
    ii = 1;
    params(ii).name = sprintf('baseline (log utility)'); 
    params(ii).OneAsset = 0;
    params(ii).DirIncomeProcess = 'input/IncomeGrids/twopoint_3_5';
    params(ii).chi0 = 0;
    params(ii).chi1 = 0.15;
    params(ii).chi2 = 0.25;
    params(ii).a_lb = 0.25;
    params(ii).riskaver = 1;
    params(ii).invies = 1;
    params(ii).SDU = 0;
    params(ii).r_a = 0.022866;
    params(ii).delta_HJB = 10;
    params(ii).delta_KFE = 1e6;
    params(ii).maxit_HJB = 1e6;
    params(ii).maxit_KFE = 1e6;
    params(ii).NoRisk = 0;
    params(ii).nb = 50;
    params(ii).nb_KFE = 50;
    params(ii).na = 50;
    params(ii).na_KFE = 50;
    params(ii).deathrate = 0;
    params(ii).rho = rho_ies1;
    params(ii).rhoL = 0.022;
    params(ii).transfer = 0.0081 * 2.0;
    params(ii).implicit = 0;
    params(ii).SaveResults = 0;

    %%--------------------------------------------------------------
    % WITH RETURNS RISK
    % --------------------------------------------------------------
    % illiquid_returns = linspace(0.03/4, 0.12/4, 10);
    ies_vals = [1, 1.5];
    risk_avers = [1, 2, 5, 10, 20];
    sdrs = [0, 0.01, 0.02, 0.05, 0.1, 0.15];
    chi1s = [0.15, 0.4096];
    
    RA5calibration = 1;

    ii = 2;
    for chi1 = chi1s
        for ies = ies_vals
            for risk_aver = risk_avers
                for sd_r = sdrs
                    params(ii).name = sprintf('SDU with riskaver%f, sigma_r%f, ies%f', risk_aver, sd_r, ies); 
                    params(ii).OneAsset = 0;
                    params(ii).DirIncomeProcess = 'input/IncomeGrids/twopoint_3_5';
                    params(ii).chi0 = 0;
                    params(ii).chi1 = chi1;
                    params(ii).chi2 = 0.25;
                    params(ii).a_lb = 0.25;
                    params(ii).riskaver = risk_aver;
                    params(ii).invies = 1 / ies;
                    params(ii).SDU = 1;
                    params(ii).maxit_HJB = 1e6;
                    params(ii).maxit_KFE = 1e6;
                    params(ii).crit_HJB = 1e-9;
                    params(ii).sigma_r = sd_r;
                    params(ii).retrisk_KFE = 0;
                    params(ii).NoRisk = 0;
                    params(ii).delta_HJB = 10;
                    params(ii).delta_KFE = 1e6;
                    params(ii).nb = 50;
                    params(ii).nb_KFE = 50;
                    params(ii).na = 50;
                    params(ii).na_KFE = 50;
                    params(ii).deathrate = 0;
                    params(ii).implicit = 0;
                    params(ii).transfer = 0.0081 * 2.0;
                    params(ii).r_b = 0.02 / 4;
                    params(ii).iterateKFE = 1;
                    params(ii).SaveResults = 0;

                    % discount factor
                    if ies == 1
                        params(ii).rho = rho_ies1;
                    else
                        params(ii).rho = rho_ies1_5;
                    end

                    % risk_aver = 1 special case not coded for ies ~= 1
                    if (ies == 1.5) && (risk_aver == 1)
                        params(ii).riskaver = 1.01;
                    end
                    
                    % log utility case
                    if (risk_aver == 1) && (ies == 1)
                        params(ii).SDU = 0;
                    end

                    if RA5calibration == 0
                        % set delta_HJB to depend on parameters
                        if ies == 1
                            if (risk_aver == 5) && (sd_r > 0.01)
                                params(ii).delta_HJB = 0.2;
                            elseif (risk_aver == 10) && (sd_r <= 0.01)
                                params(ii).delta_HJB = 1;
                            elseif risk_aver == 10
                                params(ii).delta_HJB = 0.2;
                            elseif (risk_aver == 20) && (sd_r <= 0.01)
                                params(ii).delta_HJB = 0.1;
                            elseif (risk_aver == 20) && (sd_r <= 0.1)
                                params(ii).delta_HJB = 0.025;
                            elseif (risk_aver == 20) && (sd_r < 0.15)
                                params(ii).delta_HJB = 0.01;
                            elseif (risk_aver == 20)
                                params(ii).delta_HJB = 0.005;
                            end
                        elseif ies == 1.5
                            if risk_aver <= 2
                                params(ii).delta_HJB = 2;
                            elseif (risk_aver == 5) && (sd_r <= 0.02)
                                params(ii).delta_HJB = 2;
                            elseif (risk_aver == 5)
                                params(ii).delta_HJB = 0.5;
                            elseif (risk_aver == 10) && (sd_r <= 0.02)
                                params(ii).delta_HJB = 0.5;
                            elseif risk_aver == 10
                                params(ii).delta_HJB = 0.1;
                            elseif (risk_aver == 20) && (sd_r <= 0.05)
                                params(ii).delta_HJB = 0.05;
                            elseif (risk_aver == 20) && (sd_r < 0.15)
                                params(ii).delta_HJB = 0.01;
                            elseif (risk_aver == 20)
                                params(ii).delta_HJB = 0.005;
                            end
                        end
                    else
                        if ies == 1
                            if (risk_aver == 5) && (sd_r > 0.01)
                                params(ii).delta_HJB = 0.2;
                            elseif (risk_aver == 10) && (sd_r <= 0.01)
                                params(ii).delta_HJB = 1;
                            elseif risk_aver == 10
                                params(ii).delta_HJB = 0.2;
                            elseif (risk_aver == 20) && (sd_r <= 0.01)
                                params(ii).delta_HJB = 0.1;
                            elseif (risk_aver == 20) && (sd_r <= 0.1)
                                params(ii).delta_HJB = 0.025;
                            elseif (risk_aver == 20) && (sd_r < 0.15)
                                params(ii).delta_HJB = 0.01;
                            elseif (risk_aver == 20)
                                params(ii).delta_HJB = 0.005;
                            end
                        elseif ies == 1.5
                            if risk_aver <= 2
                                params(ii).delta_HJB = 2;
                            elseif (risk_aver == 5) && (sd_r <= 0.02)
                                params(ii).delta_HJB = 2;
                            elseif (risk_aver == 5)
                                params(ii).delta_HJB = 0.5;
                            elseif (risk_aver == 10) && (sd_r <= 0.02)
                                params(ii).delta_HJB = 0.5;
                            elseif risk_aver == 10
                                params(ii).delta_HJB = 0.1;
                            elseif (risk_aver == 20) && (sd_r <= 0.05)
                                params(ii).delta_HJB = 0.05;
                            elseif (risk_aver == 20) && (sd_r < 0.15)
                                params(ii).delta_HJB = 0.01;
                            elseif (risk_aver == 20)
                                params(ii).delta_HJB = 0.005;
                            end
                        end
                    end

                    ii = ii + 1;
                end
            end
        end
    end

    %% DO NOT CHANGE BELOW

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = setup.Params(runopts,chosen_param);

end