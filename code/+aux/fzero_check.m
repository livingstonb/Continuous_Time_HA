function stop = fzero_check(~,optimValues,~,p)
    % This is an output function used with fzero to set a maximum value for
    % the number of iterations. After hitting this value, the stop flag is
    % passed as true.
    
    if optimValues.funccount >= p.maxit_AY
        stop = true;
    else
        stop = false;
    end
end