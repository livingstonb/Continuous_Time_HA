function display_exception_stack(me)
	disp(me.message)

	filestack = {me.stack.file};
	namestack = {me.stack.name};
	linestack = [me.stack.line];
	
	for istack = 1:numel(me.stack)
		fprintf("File: %s, Name: %s, Line: %d\n",...
			me.stack(istack).file, me.stack(istack).name, me.stack(istack).line)
	end
end