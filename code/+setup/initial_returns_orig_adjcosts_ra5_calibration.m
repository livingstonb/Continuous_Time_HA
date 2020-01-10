function [r_b_0, r_a_0] = initial_returns_orig_adjcosts_ra5_calibration(p)

	if p.invies == 1
	    switch p.riskaver
	        case 1
	            switch p.sigma_r
	                case 0
	                    r_b_0 = 0.022;
	                    r_a_0 = 0.061866;
	                case 0.01
	                    r_b_0 = 0.005006;
	                    r_a_0 = 0.022877;
	                case 0.02
	                    r_b_0 = 0.005026;
	                    r_a_0 = 0.022909;
	                case 0.05
	                    r_b_0 = 0.005154;
	                    r_a_0 = 0.02313;
	                case 0.1
	                    r_b_0 = 0.005539;
	                    r_a_0 = 0.023891;
	                case 0.15
	                    r_b_0 = 0.006049;
	                    r_a_0 = 0.025091;
	            end
	        case 2
	            switch p.sigma_r
	                case 0
	                    r_b_0 = 0.045451;
	                    r_a_0 = 0.059339;
	                case 0.01
	                    r_b_0 = -0.009856;
	                    r_a_0 = 0.019937;
	                case 0.02
	                    r_b_0 = -0.00982;
	                    r_a_0 = 0.020134;
	                case 0.05
	                    r_b_0 = -0.0098;
	                    r_a_0 = 0.022134;
	                case 0.1
	                    r_b_0 = -0.010709;
	                    r_a_0 = 0.024941;
	                case 0.15
	                    r_b_0 = -0.00961;
	                    r_a_0 = 0.032941;
	            end
	        case 5
	            switch p.sigma_r
	                case 0
	                    r_b_0 = -0.051442;
	                    r_a_0 = 0.009666;
	                case 0.01
	                    r_b_0 = -0.051958;
	                    r_a_0 = 0.009898;
	                case 0.02
	                    r_b_0 = -0.05323;
	                    r_a_0 = 0.010561;
	                case 0.05
	                    r_b_0 = -0.057679;
	                    r_a_0 = 0.014174;
	                case 0.1
	                    r_b_0 = -0.065535;
	                    r_a_0 = 0.022612;
	                case 0.15
	                    r_b_0 = -0.07403;
	                    r_a_0 = 0.03228;
	            end
	        case 10
	            switch p.sigma_r
	                case 0
	                    r_b_0 = -0.092879;
	                    r_a_0 = 0.00526;
	                case 0.01
	                    r_b_0 = -0.094133;
	                    r_a_0 = 0.005732;
	                case 0.02
	                    r_b_0 = -0.098433;
	                    r_a_0 = 0.009332;
	                case 0.05
	                    r_b_0 = -0.105433;
	                    r_a_0 = 0.015832;
	                case 0.1
	                    r_b_0 = -0.118208;
	                    r_a_0 = 0.02589;
	                case 0.15
	                    r_b_0 = -0.131433;
	                    r_a_0 = 0.036032;
	            end
	        case 20
	            switch p.sigma_r
	                case 0
	                    r_b_0 = -0.124384;
	                    r_a_0 = 0.004371;
	                case 0.01
	                    r_b_0 = -0.128384;
	                    r_a_0 = 0.006371;
	                case 0.02
	                    r_b_0 = -0.133384;
	                    r_a_0 = 0.012371;
	                case 0.05
	                    r_b_0 = -0.144384;
	                    r_a_0 = 0.018371;
	                case 0.1
	                    r_b_0 = -0.142;
	                    r_a_0 = 0.026471;
	                case 0.15
	                    r_b_0 = -0.1689;
	                    r_a_0 = 0.038;
	            end
	    end
	else
	    switch p.riskaver
	        case 1.01
	            switch p.sigma_r
	                case 0
	                    r_b_0 = 0.005000;
	                    r_a_0 = 0.019368;
	                case 0.01
	                    r_b_0 = 0.00503;
	                    r_a_0 = 0.019402;
	                case 0.02
	                    r_b_0 = 0.005118;
	                    r_a_0 = 0.019501;
	                case 0.05
	                    r_b_0 = 0.005154;
	                    r_a_0 = 0.0199;
	                case 0.1
	                    r_b_0 = 0.005539;
	                    r_a_0 = 0.021891;
	                case 0.15
	                    r_b_0 = 0.006049;
	                    r_a_0 = 0.023091;
	            end
	        case 2
	            switch p.sigma_r
	                case 0
	                    r_b_0 = -0.009868;
	                    r_a_0 = 0.01987;
	                case 0.01
	                    r_b_0 = -0.009856;
	                    r_a_0 = 0.019937;
	                case 0.02
	                    r_b_0 = -0.00982;
	                    r_a_0 = 0.020134;
	                case 0.05
	                    r_b_0 = -0.0098;
	                    r_a_0 = 0.022134;
	                case 0.1
	                    r_b_0 = -0.010709;
	                    r_a_0 = 0.024941;
	                case 0.15
	                    r_b_0 = -0.00961;
	                    r_a_0 = 0.032941;
	            end
	        case 5
	            switch p.sigma_r
	                case 0
	                    r_b_0 = -0.051442;
	                    r_a_0 = 0.009666;
	                case 0.01
	                    r_b_0 = -0.051958;
	                    r_a_0 = 0.009898;
	                case 0.02
	                    r_b_0 = -0.05323;
	                    r_a_0 = 0.010561;
	                case 0.05
	                    r_b_0 = -0.057679;
	                    r_a_0 = 0.014174;
	                case 0.1
	                    r_b_0 = -0.065535;
	                    r_a_0 = 0.022612;
	                case 0.15
	                    r_b_0 = -0.07403;
	                    r_a_0 = 0.03228;
	            end
	        case 10
	            switch p.sigma_r
	                case 0
	                    r_b_0 = -0.08399;
	                    r_a_0 = 0.00631;
	                case 0.01
	                    r_b_0 = -0.084389;
	                    r_a_0 = 0.006856;
	                case 0.02
	                    r_b_0 = -0.084691;
	                    r_a_0 = 0.006332;
	                case 0.05
	                    r_b_0 = -0.085239;
	                    r_a_0 = 0.008287;
	                case 0.1
	                    r_b_0 = -0.097537;
	                    r_a_0 = 0.031703;
	                case 0.15
	                    r_b_0 = -0.105433;
	                    r_a_0 = 0.041032;
	            end
	        case 20
	            switch p.sigma_r
	                case 0
	                    r_b_0 = -0.101910;
	                    r_a_0 = 0.007058;
	                case 0.01
	                    r_b_0 = -0.103299;
	                    r_a_0 = 0.007687;
	                case 0.02
	                    r_b_0 = -0.104012;
	                    r_a_0 = 0.010050;
	                case 0.05
	                    r_b_0 = -0.107618;
	                    r_a_0 = 0.028151;
	                case 0.1
	                    r_b_0 = -0.148;
	                    r_a_0 = 0.030471;
	                case 0.15
	                    r_b_0 = -0.1689;
	                    r_a_0 = 0.042;
	            end
	    end
	end
end