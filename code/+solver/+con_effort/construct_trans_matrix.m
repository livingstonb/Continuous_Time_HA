function A = construct_trans_matrix(p,income,grids,policies)

    % A or Au built without income and death transitions
    % ------------- unpack objects-----------------
    nb = numel(grids.b.vec);
    nc = numel(grids.c.vec);
    ny = income.ny;
    
    c = policies.c;
    s = policies.s;
    h = policies.h;
    % ---------------------------------------------- 
    
    if strcmp(p.hdef,"cdot/c")
        cdot = h .* c;
    elseif ismember(p.hdef,{'cdot','inf'})
        cdot = h;
    end

    %% ---------------------------------------------
    % SAVING
    % ----------------------------------------------

    sneg = min(s,0);
    sneg(1,:,:) = 0;
    spos = max(s,0);
    spos(nb,:,:) = 0;

    % lower diagonal
    xx = - sneg ./ grids.b.dB;
    xx = xx(:);
    xx = xx(2:end); % force first entry to be omitted in matrix
    xx(end+1) = 0; % needs to be consistent length

    % upper diagonal
    zz = spos ./ grids.b.dF;
    zz = zz(:);
    zz = [0;zz(1:end-1)]; % first entry automatically omitted so add a zero

    % main diagonal
    yy = - spos./grids.b.dF + sneg./grids.b.dB;
            
    yy = yy(:);

    A = spdiags([xx yy zz],[-1 0 1],nb*nc*ny,nb*nc*ny);

    %% ---------------------------------------------
    % CONSUMPTION
    % ----------------------------------------------
    cdotneg = min(cdot,0);
    cdotneg(:,1,:) = 0;
    cdotpos = max(cdot,0);
    cdotpos(:,nc,:) = 0;

    % lower diag
    clower = - cdotneg ./ grids.c.dB;
    clower = clower(:);
    clower = clower(nb+1:end);
    clower = [clower(:);zeros(nb,1)];
    
    % main diag
    cmain = - cdotpos./grids.c.dF + cdotneg./grids.c.dB;
    cmain = cmain(:);
    
    % upper diag
    cupper = cdotpos ./ grids.c.dF;
    cupper = cupper(:);
    cupper = [zeros(nb,1);cupper(1:end-nb)];
    
    A = A + spdiags([clower cmain cupper],[-nb 0 nb],nb*nc*ny,nb*nc*ny);