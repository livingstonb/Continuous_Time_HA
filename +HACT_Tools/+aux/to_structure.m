function s = to_structure(objs)
    % Converts an object to a structure
    % by copying its fields. Not guaranteed
    % to work on all objects.
    
    ofields = fields(objs);
    for is = 1:numel(objs)
        for ifield = 1:numel(ofields)
            s(is).(ofields{ifield}) = objs(is).(ofields{ifield});
        end
    end
end