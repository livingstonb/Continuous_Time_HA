function outparams = chi0_tests(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,...

    chi0s = [0.001,0.0025,0.005];
    chi1s = [0.1,0.2,0.5];
    
    ras_chi0_001 = [0.0282004169674092;0.0253830746362491;0.0203726112664956];
        
    ras_chi0_0025 = [0.0282254647550881;0.0254205136551605;0.0204459253227402];
    
    ras_chi0_005 = [0.0282668949847884;0.0254982338526208;0.0205676145496176];
    
    ii = 1;
	for chi0 = chi0s
    for chi1 = chi1s
    	chi2 = 0.5;
        
	    params(ii).name = sprintf('contB, chi0_%0.2f',chi0); 
	    params(ii).OneAsset = 0;
	    params(ii).DirIncomeProcess = 'IncomeGrids/continuous_b';
	    params(ii).chi2 = chi2;
	    params(ii).chi1 = chi1;
	    params(ii).chi0 = chi0;
	    params(ii).a_lb = 0.25;
        
        j = mod(ii,3);
        if j == 0
            j = 3;
        end
        
        switch chi0
            case 0.001
                params(ii).r_a = ras_chi0_001(j);
            case 0.0025
                params(ii).r_a = ras_chi0_0025(j);
            case 0.005
                params(ii).r_a = ras_chi0_005(j);
        end
       
        ii = ii + 1;
    end
    end

    
 %    %% continuous_b, chi2 = 0.5
    
 %    ras = [0.0315337763625455
 %        0.0301320288894038
 %        0.0293393022213278
 %        0.0287588072381727
 %        0.0282538595226361];
        
 %    for eps = [0]
	% for ip = 1:numel(chi0s)
 %    	chi1 = 0.1;
 %    	chi2 = 0.5;
 %    	ra = ras(ip) + eps;
        
	%     params(end+1).name = sprintf('contB, chi0_%0.2f',chi0); 
	%     params(end).OneAsset = 0;
	%     params(end).DirIncomeProcess = 'IncomeGrids/continuous_b';
	%     params(end).chi2 = chi2;
	%     params(end).chi1 = chi1;
	%     params(end).chi0 = 0;
	%     params(end).a_lb = 0.25;
	%     params(end).r_a = ra;
 %    end
 %    end
    
 %    %% 0.25
 %    chi2 = 0.25;
 %    ras = [0.0276395376993103;0.0272237189969445;0.0268694847830224;0.0265623593325322;0.0262878669204621;0.0260339044482922;0.0257978593234693;0.0255787409161354;0.0253725559412461;0.0251760086425505;0.0249915960341563;0.0248224553905225;0.0246684922571530;0.0245269669546823;0.0243877906332216;0.0242426899055111;0.0240980024053304;0.0239637185217187;0.0238498286437152];

 %    chi1s = 0.1:0.05:1;
 %    for add = [0]
 %        for i = 1:numel(chi1s)

 %            ra = ras(i) + add;
 %            chi1 = chi1s(i);

 %            params(end+1).name = ['chi1_',num2str(chi1),', chi2_',num2str(chi2)];
 %            params(end).chi0 = 0;
 %            params(end).chi1 = chi1;
 %            params(end).chi2 = chi2;
 %            params(end).a_lb = 0.25;
 %            params(end).r_a = ra;
 %        end
 %    end

 %    %% Chi2 = 0.333
 %    chi2 = 0.333;
 %    ras = [0.0266352878553511;0.0260286877288234;0.0254769652937515;0.0249753112792640;0.0245189164144894;0.0241031000752887;0.0237236962244517;0.0233766816546957;0.0230580898915175;0.0227641823929926;0.0224920756147316;0.0222391313300261;0.0220028375853561;0.0217811395355858;0.0215736844959291;0.0213796395979645;0.0211945490783803;0.0210130514501425;0.0208297852262168];
 %    chi1s = 0.1:0.05:1;
    
    
 %    for add = [0]
 %        for i = 1:numel(chi1s)

 %            ra = ras(i) + add;
 %            chi1 = chi1s(i);

 %            params(end+1).name = ['chi1_',num2str(chi1),', chi2_',num2str(chi2)];
 %            params(end).chi0 = 0;
 %            params(end).chi1 = chi1;
 %            params(end).chi2 = chi2;
 %            params(end).a_lb = 0.25;
 %            params(end).r_a = ra;
 %        end
 %    end
    
 %    %% Chi2 = 1
 %    chi2 = 1;
 %    ras = [0.0194518644146336;0.0167851960583025;0.0148971673720334;0.0135828543860802;0.0126373331306967;0.0118854700490672;0.0112712932360973;0.0107685521641382;0.0103507201636011;0.00999406159380381;0.00968628107163059;0.00941841070482711;0.00918335230701831;0.00897495710612281;0.00878900428135457;0.00862199053754181;0.00847135473067441;0.00833477125453265;0.00820991450289681];
 %    chi1s = 0.1:0.05:1;
 %    for add = [0]
 %        for i = 1:numel(chi1s)        
 %            ra = ras(i) + add;
 %            chi1 = chi1s(i);

 %            params(end+1).name = ['chi1_',num2str(chi1),', chi2_',num2str(chi2)];
 %            params(end).chi0 = 0;
 %            params(end).chi1 = chi1;
 %            params(end).chi2 = chi2;
 %            params(end).a_lb = 0.25;
 %            params(end).r_a = ra;
 %        end
 %    end

 %    %% Chi2 = 0.5
	% chi2 = 0.5;
 %    ras = [0.0248644665224153;0.0235921222145708;0.0224957838107575;0.0215560181290627;0.0207533919875736;0.0200676587350894;0.0194753178432566;0.0189543192612625;0.0184916687256102;0.0180772444246925;0.0177033585671444;0.0173633326409063;0.0170520912308980;0.0167654596136322;0.0165012627350145;0.0162569927926730;0.0160268113217330;0.0158040471916944;0.0155820292720568];
 %    chi1s = 0.1:0.05:1;
 %    for add = [0]
 %        for i = 1:numel(chi1s)
 %            ra = ras(i) + add;
 %            chi1 = chi1s(i);

 %            params(end+1).name = ['chi1_',num2str(chi1),', chi2_',num2str(chi2)];
 %            params(end).chi0 = 0;
 %            params(end).chi1 = chi1;
 %            params(end).chi2 = chi2;
 %            params(end).a_lb = 0.25;
 %            params(end).r_a = ra;
 %        end
 %    end

 %    %% Chi2 = 2
 %    chi2 = 2;
 %    ras = [0.0133868407289295;0.0108261398939028;0.00918483987201127;0.00819889855944732;0.00760427385240319;0.00717758335245619;0.00685808348272363;0.00662458169844137;0.00645145076177941;0.00631535887964103;0.00620659073092848;0.00611909324380778;0.00604784587150111;0.00598853074551334;0.00593860818542389;0.00589615204425937;0.00585991212076053;0.00582880720009666;0.00580175606743705];
 %    chi1s = 0.1:0.05:1;
 %    for add = [0]
 %        for i = 1:numel(chi1s)
 %            ra = ras(i) + add;
 %            chi1 = chi1s(i);

 %            params(end+1).name = ['chi1_',num2str(chi1),', chi2_',num2str(chi2)];
 %            params(end).chi0 = 0;
 %            params(end).chi1 = chi1;
 %            params(end).chi2 = chi2;
 %            params(end).a_lb = 0.25;
 %            params(end).r_a = ra;
 %        end
 %    end

    %%

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = model_objects.Params(runopts,chosen_param);

end
