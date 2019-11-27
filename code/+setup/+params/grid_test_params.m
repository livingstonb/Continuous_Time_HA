function outparams = grid_test_params(runopts)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,...
    params(1).name = 'n = 100, nKFE = 100'; 
    
%     params(2).name = 'low amax';
%     params(2).amax = 100;
%     
%     params(3).name = 'high amax';
%     params(3).amax = 1000;
%     
%     params(4).name = 'low bmax';
%     params(4).bmax = 50;
%     
%     params(5).name = 'high bmax';
%     params(5).bmax = 200;
    
    params(end+1).name = 'low bmax and amax';
    params(end).bmax = 50;
    params(end).amax = 100;
    
    params(end+1).name = 'high bmax and amax';
    params(end).bmax = 200;
    params(end).amax = 1000;
    
    params(end+1).name = 'low grid curv';
    params(end).b_gcurv_pos = 0.2;
    params(end).b_gcurv_neg = 0.2;
    params(end).a_gcurv = 0.2;
    
    params(end+1).name = 'high nb and na';
    params(end).nb = 150;
    params(end).na = 150;
    params(end).nb_KFE = 120;
    params(end).na_KFE = 120;
% 
%     params(2).name = 'n = 100, nKFE = 200';
%     params(2).nb_KFE = 200;
%     params(2).na_KFE = 200;
% 
%     params(3).name = 'n = 200, nKFE = 100';
%     params(3).nb = 200;
%     params(3).na = 200;
% 
%     params(4).name = 'n = 150, nKFE = 150'
%     params(4).nb = 150;
%     params(4).nb_KFE = 150;
%     params(4).na = 150;
%     params(4).na_KFE = 150;
% 
%     params(5).name = 'n = 200, nKFE = 200';
%     params(5).nb = 200;
%     params(5).na = 200;
%     params(5).nb_KFE = 200;
%     params(5).na_KFE = 200;

    % ---------------------------------------------------------------------

    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACTLib.model_objects.Params(runopts,chosen_param);

end
