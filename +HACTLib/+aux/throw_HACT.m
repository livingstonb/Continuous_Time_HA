function create_exception
	

end

function msgId = gen_msgId(caller, varargin)
	msgId = "HACTLib";
	for k = 1:nargin-1
		msgId = cat(msgId, ":", varargin{k});
	end
end