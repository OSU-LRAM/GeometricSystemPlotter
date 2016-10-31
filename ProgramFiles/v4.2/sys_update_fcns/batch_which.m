function files = batch_which(target)
% Batch application of the which command to a cell array of filename
% strings

	if ischar(target)
		
		files = pathlookup(target);
		
	elseif iscellstr(target)

		files = cellfun(@(t) pathlookup(t),target,'UniformOutput',false);
	
	else
		
		error('batch_which argument must be a string or cell array of strings')
		
	end
	
end


function pathout = pathlookup(pathin)

		if exist(pathin,'dir')
			
			pathout = pathin;
			
		else
		
			pathout = which(pathin);
			
		end
		
end