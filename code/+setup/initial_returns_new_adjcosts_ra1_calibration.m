function [r_b_0, r_a_0] = initial_returns_new_adjcosts_ra1_calibration(p)

	if p.invies == 1
	    switch p.riskaver
	        case 1
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = 0.005;
	                    r_a_0 = 0.022612;
	                case 0.01
	                    converged = true;
	                    r_b_0 = 0.005006;
	                    r_a_0 = 0.022627;
	                case 0.02
	                    converged = true;
	                    r_b_0 = 0.005026;
	                    r_a_0 = 0.020909;
	                case 0.05
	                    converged = true;
	                    r_b_0 = 0.005054;
	                    r_a_0 = 0.022669;
	                case 0.1
	                    converged = true;
	                    r_b_0 = 0.005339;
	                    r_a_0 = 0.023991;
	                case 0.15
	                    converged = true;
	                    r_b_0 = 0.005408;
	                    r_a_0 = 0.025668;
	            end
	        case 2
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = -0.006883;
	                    r_a_0 = 0.019505;
	                case 0.01
	                    converged = true;
	                    r_b_0 = -0.006915;
	                    r_a_0 = 0.019600;
	                case 0.02
	                    converged = true;
	                    r_b_0 = -0.007010;
	                    r_a_0 = 0.019881;
	                case 0.05
	                    converged = true;
	                    r_b_0 = -0.007469;
	                    r_a_0 = 0.021652;
	                case 0.1
	                    converged = true;
	                    r_b_0 = -0.009050;
	                    r_a_0 = 0.026614;
	                case 0.15
	                    converged = true;
	                    r_b_0 = -0.013651; % converged, error with no risk
	                    r_a_0 = 0.033504; % converged, error with no risk
	            end
	        case 5
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = -0.040727;
	                    r_a_0 = 0.008936;
	                case 0.01
	                    converged = true;
	                    r_b_0 = -0.041411;
	                    r_a_0 = 0.009220;
	                case 0.02
	                    r_b_0 = -0.043032;
	                    r_a_0 = 0.011039;
	                case 0.05
	                    r_b_0 = -0.047679;
	                    r_a_0 = 0.015174;
	                case 0.1
	                    r_b_0 = -0.053535;
	                    r_a_0 = 0.021612;
	                case 0.15
	                    r_b_0 = -0.063892;
	                    r_a_0 = 0.030143;
	            end
	        case 10
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = -0.076174;
	                    r_a_0 = 0.003779;
	                case 0.01
	                    converged = true;
	                    r_b_0 = -0.078166;
	                    r_a_0 = 0.004316;
	                case 0.02
	                    r_b_0 = -0.082166;
	                    r_a_0 = 0.005116;
	                case 0.05
	                    r_b_0 = -0.090433;
	                    r_a_0 = 0.006132;
	                case 0.1
	                    r_b_0 = -0.099433;
	                    r_a_0 = 0.007032;
	                case 0.15
	                    r_b_0 = -0.111433;
	                    r_a_0 = 0.008032;
	            end
	        case 20
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = -0.104539;
	                    r_a_0 = 0.002535;
	                case 0.01
	                    r_b_0 = -0.107539;
	                    r_a_0 = 0.003371;
	                case 0.02
	                    r_b_0 = -0.111539;
	                    r_a_0 = 0.004071;
	                case 0.05
	                    r_b_0 = -0.114384;
	                    r_a_0 = 0.005371;
	                case 0.1
	                    r_b_0 = -0.1252;
	                    r_a_0 = 0.006971;
	                case 0.15
	                    r_b_0 = -0.1289;
	                    r_a_0 = 0.012;
	            end
	    end
	else
	    switch p.riskaver
	        case 1.01
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = 0.005;
	                    r_a_0 = 0.019791;
	                case 0.01
	                    converged = true;
	                    r_b_0 = 0.005026;
	                    r_a_0 = 0.019836;
	                case 0.02
	                    converged = true;
	                    r_b_0 = 0.005104;
	                    r_a_0 = 0.019966;
	                case 0.05
	                    converged = true;
	                    r_b_0 = 0.005612;
	                    r_a_0 = 0.020795;
	                case 0.1
	                    r_b_0 = 0.0063308;
	                    r_a_0 = 0.023997;
	                case 0.15
	                    r_b_0 = 0.0071308;
	                    r_a_0 = 0.030597;
	            end
	        case 2
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = -0.005964;
	                    r_a_0 = 0.017394;
	                case 0.01
	                    converged = true;
	                    r_b_0 = -0.005943;
	                    r_a_0 = 0.017522;
	                case 0.02
	                    converged = true;
	                    r_b_0 = -0.005874;
	                    r_a_0 = 0.017893;
	                case 0.05
	                    converged = true;
	                    r_b_0 = -0.005546;
	                    r_a_0 = 0.020098;
	                case 0.1
	                    r_b_0 = -0.007046;
	                    r_a_0 = 0.0281971;
	                case 0.15
	                    converged = true;
	                    r_b_0 = -0.009953;
	                    r_a_0 = 0.034458;
	            end
	        case 5
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = -0.038573;
	                    r_a_0 = 0.008789;
	                case 0.01
	                    converged = true;
	                    r_b_0 = -0.039105;
	                    r_a_0 = 0.009146;
	                case 0.02
	                    converged = true;
	                    r_b_0 = -0.040288;
	                    r_a_0 = 0.010153;
	                case 0.05
	                    r_b_0 = -0.043872;
	                    r_a_0 = 0.020275;
	                case 0.1
	                    r_b_0 = -0.051544;
	                    r_a_0 = 0.030421;
	                case 0.15
	                    converged = true;
	                    r_b_0 = -0.061004;
	                    r_a_0 = 0.045998;
	            end
	        case 10
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = -0.070741;
	                    r_a_0 = 0.005279;
	                case 0.01
	                    converged = true;
	                    r_b_0 = -0.071900;
	                    r_a_0 = 0.005972;
	                case 0.02
	                    r_b_0 = -0.074081;
	                    r_a_0 = 0.010808;
	                case 0.05
	                    converged = true;
	                    r_b_0 = -0.082419;
	                    r_a_0 = 0.02236;
	                case 0.1
	                    r_b_0 = -0.092600;
	                    r_a_0 = 0.039868;
	                case 0.15
	                    converged = true;
	                    r_b_0 = -0.1099149;
	                    r_a_0 = 0.059595;
	            end
	        case 20
	            switch p.sigma_r
	                case 0
	                    converged = true;
	                    r_b_0 = -0.092393;
	                    r_a_0 = 0.004989;
	                case 0.01
	                    r_b_0 = -0.094769;
	                    r_a_0 = 0.005029;
	                case 0.02
	                    r_b_0 = -0.100461;
	                    r_a_0 = 0.006972;
	                case 0.05
	                    r_b_0 = -0.110461;
	                    r_a_0 = 0.011272;
	                case 0.1
	                    r_b_0 = -0.114661;
	                    r_a_0 = 0.04572;
	                case 0.15
	                    r_b_0 = -0.127261;
	                    r_a_0 = 0.075572;
	            end
	    end
	end
end