function cell_out = get_all_values(struct_array,field1,field1Index,field2,field2Index,field3,field3Index)
    % produces cell array of the form:
    % {struct_array(1).field1.field2(field2Index,...,struct_array(n).field1.field2(field2Index)}
    %
    % field1, etc... must be strings
    %
    % field1Index is an integer index for field 1

    cell_out = cell(1,numel(struct_array));

    % find number of fields to access
    numEntries = nargin - 1;
    numFields = ceil(numEntries/2);

    switch nargin
    	case 2
    		field1Index = 0;
    	case 4
    		field2Index = 0;
    	case 6
    		field3Index = 0;
	end

    switch numFields
	    case 1
	    	fields = {field1};
	    	indices = [field1Index];
	    case 2
	    	fields = {field1,field2};
	    	indices = [field1Index,field2Index];
	    case 3
	    	fields = {field1,field2,field3};
	    	indices = [field1Index,field2Index,field3Index];
	end

    numStructs = numel(struct_array);
    
    for istruct = 1:numStructs
    	cell_out{istruct} = struct_array(istruct);
    end

    for ifield = 1:numFields
    	cell_out = extract_structure(cell_out,fields{ifield},indices(ifield));
    end
end

function cellArrayOutput = extract_structure(...
	cellArrayOfStructures,field,fieldIndex)
	
	if ~exist('fieldIndex','var')
		fieldIndex = 0;
	end

	n = numel(cellArrayOfStructures);
	cellArrayOutput = cell(1,n);

	for ii = 1:n
		ii_array = cellArrayOfStructures{ii}.(field);
        if fieldIndex == 0
            cellArrayOutput{ii} = ii_array;
        else
            cellArrayOutput{ii} = ii_array(fieldIndex);
        end
	end
end