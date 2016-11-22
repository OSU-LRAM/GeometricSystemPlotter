function unmet_dependencies = depcheck(target,varargin)
% Check if the target function is older than any file it depends on, or
% younger than any file it produces. varargin should contain any other
% arguments passed to the target file to help it in its dependency-check


	% First, if the target function is an m-file, query it for any
	% dependencies
	if strcmp('.m',target(end-1:end))
		
		try 
			
			[depstructure,full_target] = parse_dependencies(target,varargin{:});

		catch
			depstructure = [];
			full_target = fullfile(pwd,target);
		end
		
	else % file is not a function, return that we don't know any unmet dependencies for it
		
		unmet_dependencies = 0;
		return
		
	end
	
	
	% Check if any of the dependencies or products identified for the
	% target are of an age to require that the target be rerun
	if isfield(depstructure,'product') && ~isempty(depstructure.product)
		
		unmet_dependencies_1 = check_if_older(depstructure.product,full_target);
				
		if isfield(depstructure,'dependency') && ~isempty(depstructure.dependency)

			unmet_dependencies_2 = check_if_older(depstructure.product,depstructure.dependency);
			
		else
			
			unmet_dependencies_2 = 0;
			
		end
		
		% check if any dependency conditions were met
		unmet_dependencies = unmet_dependencies_1 || unmet_dependencies_2;	

	else
	% If the target is not enabled for depcheck, or returns no
	% dependencies or products then return that it is involved in no
	% (known) unmet dependencies
		
		unmet_dependencies = 0;
		
	end
	

	
	% If there were no first-order unmet dependencies, check the
	% dependencies of each file current target is dependent on
	if unmet_dependencies
		
		return
		
	else
		
		if isfield(depstructure,'dependency')
		
			for i = 1:length(depstructure.dependency)
				
				unmet_dependencies = unmet_dependencies || depcheck(depstructure.dependency{i});
				
			end
			
		end
		
	end
	
end