function D = batch_dir(target)
% Get the directory information on a cell array of filenames

	if ischar(target)
		
		D = dir(target);
		
		% Remove the first two entries (for . and .. ) from directory lists
		if isdir(target)
			D = D(3:end);
		end
	
		% Drill down along any directories
		subD = [];
		for i = 1:length(D)
			
			if D(i).isdir
				
				subD = cat(1,subD,batch_dir([target '/' D(i).name]));
				
			end
			
		end
		
		D = cat(1,D,subD);
	
	elseif iscellstr(target)
		
		% operate on each file in the list
		D = cellfun(@(t) batch_dir(t),target,'UniformOutput',false);
				
		% concatenate information received
		D = cat(1,D{~cellfun(@(x) isempty(x),D)});
		
	else
		
		error('batch_dir argument must be a string or cell array of strings')
		
	end


end