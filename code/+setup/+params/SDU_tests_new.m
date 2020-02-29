function outparams = SDU_tests(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,.

    import HACTLib.aux.set_shared_fields

    % rho for calibration based on baseline: 0.033940
    % rho for calibration based on riskaver = 5, sigma_r = 0.1: 0.114248

    %% calibrated rho's
    % adj cost original, calibrated to RA = 1
    rho_ies1_chi1_015 = 0.0317237836;
    % rho_ies1_5_chi1_015 = 0.026993; % not correct

    % % new adj cost, calibrated to IES = 1
    % rho_ies1_chi1_high = 0.031219; % (r_a = 0.22612)
    % rho_ies1_5_chi1_high = 0.025568;

%     % adj cost original, calibrated to RA = 5, sigma_r = 0.1
%     rho_ies1_chi1_015 = 0.114310;
%     rho_ies1_5_chi1_015 = 0.095630;

% %     adj cost new, calibrated to RA = 5, sigma_r = 0.1
%     rho_ies1_chi1_high = 0.108251;
%     rho_ies1_5_chi1_high = 0.090342;

    shared_params.nb = 40;
    shared_params.nb_KFE = 40;
    shared_params.na = 40;
    shared_params.na_KFE = 40;
    shared_params.min_grid_spacing = -1e10;
    % shared_params.mpc_shocks_dollars = [-1, -500, -5000, 1, 500, 5000];
    % shared_params.mpc_shocks = shared_params.mpc_shocks_dollars ./ 72000;
    shared_params.bmax = 20;
    shared_params.amax = 50;
    shared_params.b_gcurv_pos = 0.3;
    shared_params.a_gcurv = 0.3;
    shared_params.OneAsset = 0;
    shared_params.income_dir = 'twopoint_3_5';

    %%--------------------------------------------------------------
    % BASELINE
    % --------------------------------------------------------------
    ii = 1;
    params(ii).name = sprintf('baseline'); 
    params(ii).OneAsset = 0;
    params(ii).income_dir = 'twopoint_3_5';
    params(ii).chi0 = 0;
    params(ii).chi1 = 0.15;
    params(ii).chi2 = 0.25;
    params(ii).a_lb = 0.25;
    params(ii).riskaver = 1;
    params(ii).invies = 1;
    params(ii).SDU = false;
    params(ii).r_a = 0.02;
    params(ii).HJB_delta = 100;
    params(ii).KFE_delta = 1e6;
    params(ii).HJB_maxiters = 1e6;
    params(ii).KFE_maxiters = 1e6;
    params(ii).NoRisk = 0;
    params(ii).deathrate = 0;
    params(ii).rho = rho_ies1_chi1_015;
    params(ii).transfer = 0.0081 * 2.0;
    params(ii).HJB_implicit = false;
    params(ii).SaveResults = 1;

    %%--------------------------------------------------------------
    % WITH RETURNS RISK
    % --------------------------------------------------------------
    % illiquid_returns = linspace(0.03/4, 0.12/4, 10);
    ies_vals = [1]; % 1.5
    risk_avers = [1, 2, 5, 10, 20];
    sdrs = [0, 0.01, 0.02, 0.05, 0.1, 0.15];
    chi1s = [0.15]; % 0.4096
    
    RA5calibration = 1;

    ii = 2;
    for chi1 = chi1s
        for ies = ies_vals
            for risk_aver = risk_avers
                for sd_r = sdrs
                    params(ii).name = sprintf(...
                        'SDU with riskaver%f, sigma_r%f, ies%f, chi1_%f',...
                        risk_aver, sd_r, ies, chi1);
                    params(ii).OneAsset = 0;
                    params(ii).income_dir = 'twopoint_3_5';
                    params(ii).chi0 = 0;
                    params(ii).chi1 = chi1;
                    params(ii).chi2 = 0.25;
                    params(ii).a_lb = 0.25;
                    params(ii).riskaver = risk_aver;
                    params(ii).invies = 1 / ies;
                    params(ii).SDU = true;
                    params(ii).HJB_maxiters = 1e6;
                    params(ii).KFE_maxiters = 1e6;
                    params(ii).HJB_tol = 1e-9;
                    params(ii).sigma_r = sd_r;
                    params(ii).retrisk_KFE = 0;
                    params(ii).NoRisk = 0;
                    params(ii).HJB_delta = 100;
                    params(ii).KFE_delta = 1e6;
                    params(ii).deathrate = 0;
                    params(ii).HJB_implicit = false;
                    params(ii).transfer = 0.0081 * 2.0;
                    params(ii).KFE_iterative = true;
                    params(ii).SaveResults = 1;
                    params(ii).r_a = 0.02;

                    params(ii).rho = rho_ies1_chi1_015;

                    if chi1 == 0.15
                        switch risk_aver
                            case 1
                                params(ii).r_b = 0.005;
                            case 2
                                params(ii).r_b = -0.01;
                            case 5
                                params(ii).r_b = -0.02;
                            case 10
                                params(ii).r_b = -0.05;
                            case 20
                                params(ii).r_b = -0.1;
                        end
                    end



                    % discount factor
                    % if ies == 1 && (chi1 == 0.15)
                    %     params(ii).rho = rho_ies1_chi1_015;
                    % elseif chi1 == 0.15
                    %     params(ii).rho = rho_ies1_5_chi1_015;
                    % elseif ies == 1
                    %     params(ii).rho = rho_ies1_chi1_high;
                    % else
                    %     params(ii).rho = rho_ies1_5_chi1_high;
                    % end

                    % risk_aver = 1 special case not coded for ies ~= 1
                    if (ies == 1.5) && (risk_aver == 1)
                        params(ii).riskaver = 1.01;
                    end
                    
                    % log utility case
                    if (risk_aver == 1) && (ies == 1)
                        params(ii).SDU = false;
                    end

                    % if RA5calibration == 0
                    %     % set HJB_delta to depend on parameters
                    %     if ies == 1
                    %         if (risk_aver == 5) && (sd_r > 0.01)
                    %             params(ii).HJB_delta = 0.2;
                    %         elseif (risk_aver == 10) && (sd_r <= 0.01)
                    %             params(ii).HJB_delta = 0.2;
                    %         elseif risk_aver == 10
                    %             params(ii).HJB_delta = 0.05;
                    %         elseif (risk_aver == 20) && (sd_r <= 0.01)
                    %             params(ii).HJB_delta = 0.1;
                    %         elseif (risk_aver == 20) && (sd_r <= 0.1)
                    %             params(ii).HJB_delta = 0.025;
                    %         elseif (risk_aver == 20) && (sd_r < 0.15)
                    %             params(ii).HJB_delta = 0.01;
                    %         elseif (risk_aver == 20)
                    %             params(ii).HJB_delta = 0.005;
                    %         end
                    %     elseif ies == 1.5
                    %         if risk_aver <= 2
                    %             params(ii).HJB_delta = 2;
                    %         elseif (risk_aver == 5) && (sd_r <= 0.02)
                    %             params(ii).HJB_delta = 2;
                    %         elseif (risk_aver == 5)
                    %             params(ii).HJB_delta = 0.5;
                    %         elseif (risk_aver == 10) && (sd_r <= 0.02)
                    %             params(ii).HJB_delta = 0.5;
                    %         elseif risk_aver == 10
                    %             params(ii).HJB_delta = 0.1;
                    %         elseif (risk_aver == 20) && (sd_r <= 0.05)
                    %             params(ii).HJB_delta = 0.05;
                    %         elseif (risk_aver == 20) && (sd_r < 0.15)
                    %             params(ii).HJB_delta = 0.01;
                    %         elseif (risk_aver == 20)
                    %             params(ii).HJB_delta = 0.005;
                    %         end
                    %     end
                    % else
                    %     if ies == 1
                    %         if (risk_aver == 5) && (sd_r > 0.01)
                    %             params(ii).HJB_delta = 0.1;
                    %         elseif (risk_aver == 10) && (sd_r <= 0.01)
                    %             params(ii).HJB_delta = 1;
                    %         elseif risk_aver == 10
                    %             params(ii).HJB_delta = 0.1;
                    %         elseif (risk_aver == 20) && (sd_r <= 0.01)
                    %             params(ii).HJB_delta = 0.1;
                    %         elseif (risk_aver == 20) && (sd_r <= 0.1)
                    %             params(ii).HJB_delta = 0.025;
                    %         elseif (risk_aver == 20) && (sd_r < 0.15)
                    %             params(ii).HJB_delta = 0.01;
                    %         elseif (risk_aver == 20)
                    %             params(ii).HJB_delta = 0.005;
                    %         end
                    %     elseif ies == 1.5
                    %         if risk_aver <= 2
                    %             params(ii).HJB_delta = 2;
                    %         elseif (risk_aver == 5) && (sd_r <= 0.02)
                    %             params(ii).HJB_delta = 2;
                    %         elseif (risk_aver == 5)
                    %             params(ii).HJB_delta = 0.5;
                    %         elseif (risk_aver == 10) && (sd_r <= 0.02)
                    %             params(ii).HJB_delta = 0.5;
                    %         elseif risk_aver == 10
                    %             params(ii).HJB_delta = 0.1;
                    %         elseif (risk_aver == 20) && (sd_r <= 0.05)
                    %             params(ii).HJB_delta = 0.05;
                    %         elseif (risk_aver == 20) && (sd_r < 0.15)
                    %             params(ii).HJB_delta = 0.01;
                    %         elseif (risk_aver == 20)
                    %             params(ii).HJB_delta = 0.005;
                    %         end
                    %     end
                    % end

                    ii = ii + 1;
                end
            end
        end
    end

    params = set_shared_fields(params, shared_params);

    %% DO NOT CHANGE BELOW

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACTLib.model_objects.Params(runopts,chosen_param);

    %% ATTACH CALIBRATOR
    if runopts.calibrate
        if strcmp(outparams.name, 'baseline')
            calibrator = rho_calibrator(outparams, runopts);
        else
            calibrator = returns_calibrator(outparams, runopts);
        end

        calibrator.set_handle(outparams);
        outparams.set("calibrator", calibrator, true);
    end
end

function calibrator = rho_calibrator(p, runopts)
    import HACTLib.model_objects.AltCalibrator

    param_name = {'rho', 'r_a'};
    stat_name = {'totw', 'liqw'};
    stat_target = [3.5, 0.5];

    rho_bounds = [0.001, 0.05];
    ra_bounds = [p.r_b+1e-3, 0.1];

    calibrator = AltCalibrator(p, runopts, param_name,...
        stat_name, stat_target);

    calibrator.set_param_bounds(1, rho_bounds);
    calibrator.set_param_bounds(2, ra_bounds);
end

function calibrator = returns_calibrator(p, runopts)
    import HACTLib.model_objects.AltCalibrator

    param_name = {'r_b', 'r_a'};
    stat_name = {'totw', 'liqw'};
    stat_target = [3.5, 0.5];

    rb_bounds = [-0.2, 0.02];
    ra_bounds = [0.001, 0.1];

    calibrator = AltCalibrator(p, runopts, param_name,...
        stat_name, stat_target);

    calibrator.set_param_bounds(1, rb_bounds);
    calibrator.set_param_bounds(2, ra_bounds);
end