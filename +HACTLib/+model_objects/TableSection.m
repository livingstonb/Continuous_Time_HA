classdef TableSection
	properties (SetAccess = private)
		header;
		labels = {};
	end

	methods
		function obj = TableSection(labels, header)
			obj.labels = labels;
			obj.header = header;
		end
	end
end