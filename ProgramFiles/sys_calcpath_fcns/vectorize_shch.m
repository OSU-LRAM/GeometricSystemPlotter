%Vectorize the shape change equation
function p = vectorize_shch(p)
	
	pathspec = {'phi_def',;
				'dphi_def';
				'time_def';
				'phi_marker';
				'phi_arrows'
				'phi_res'
				'cBVI_method'};
				
	
	default_values = {@(t)ones(size(t));		% Placeholder
					@(t)ones(size(t));			% Placeholder
					[0 2*pi];	% default timespan
					[];			% no marker
					0;			% no arrows
					30;         % 30 points along the path
					'none'};	% default to no area integration for cBVI
	%%%%%%%%%
	% First, make sure that there is something in each of the pathspec
	% fields
	for i = 1:length(pathspec)
		
		% Make sure that phi_def is the first element in pathspec,
		% otherwise the handling of the other elements could break
		if i == 1 && ~strcmp(pathspec{i},'phi_def')
			error('Please make sure that phi_def is the first element of the pathspec list at the top of this file')
		end
		
		switch pathspec{i}

			
			case 'phi_def'
				
				% check for absence emptiness of field
				if ~isfield(p,pathspec{i}) || isempty(p.(pathspec{i}))

					% error if no shape change is actually specified
					error('Shape change has no function specified for phi_def. Use the No Shapechange menu option to generate a plot with no shape change overlaid')
				
				else % make sure phi_def is a double-nested cell function

					p.phi_def = doublewrap(p.phi_def,class(default_values{i}));
					
				end
					
			% If shape change derivative is not defined, use a
			% numerical approximation of it
			case 'dphi_def'

				% check for absence emptiness of field
				if ~isfield(p,pathspec{i}) || isempty(p.(pathspec{i}))
					
					% Create a double-wrapped set of numerical derivatives
					% If creating all the derivatives, do not warn user
					% about them
					for j = 1:numel(p.phi_def)
						for k = 1:numel(p.phi_def{j})
							p.dphi_def{j}{k} = @(T) jacobianest(@(TT) p.phi_def{j}{k}(TT),T);
						end
					end
					
				else
					
					% Fill in any unspecified derivatives with numerical
					% derivatives
                    
       				% Ensure double-wrapping of field
    				p.(pathspec{i}) = doublewrap(p.(pathspec{i}),class(default_values{i}));

					
					% first, strokes and sections that were only partly
					% addressed					
					for j = 1:min(numel(p.dphi_def),numel(p.phi_def)) % first fill in partially-defined derivatives

						if numel(p.dphi_def) < numel(p.phi_def)
							warning('Some (but not all) shape changes *segments* have analytically specified derivatives. The remainder have been filled in numerically; fully specify the dphi_def analytically to improve performance, or specify an a numerical derivative function in the shape change file to suppress this warning')
						end
			
						for k = (numel(p.dphi_def{j})+1):numel(p.phi_def{j})
							
							
							p.dphi_def{j}{k} = @(T) jacobianest(@(TT) p.phi_def{j}{k}(TT),T);
							
						end
						
					end
					
					% second, strokes that were completely unspecified
					if numel(p.dphi_def) < numel(p.phi_def)
						warning('Some (but not all) shape changes have analytically specified derivatives. The remainder have been filled in numerically; fully specify the dphi_def analytically to improve performance, or specify an a numerical derivative function in the shape change file to suppress this warning')
					end
					
					for j = (numel(p.dphi_def)+1) : numel(p.phi_def) % first fill in partially-defined derivatives
			
						for k = 1:numel(p.phi_def{j})
							
							p.dphi_def{j}{k} = @(T) jacobianest(@(TT) p.phi_def{j}{k}(TT),T);
							
						end
						
					end
					

					

				end
		
			otherwise
			
				% check for absence emptiness of field
				if ~isfield(p,pathspec{i}) || isempty(p.(pathspec{i}))
					
					% Create a single instance of the default object
					p.(pathspec{i}) = default_values{i};
                    
                end
					
                    % Ensure double-wrapping of field
                    p.(pathspec{i}) = doublewrap(p.(pathspec{i}),class(default_values{i}));

				
				
				%%%%
				% Match the dimension of the field to that of the phi_def
				% function
				
				% number of paths here and in phi_def
				n_thisfield = numel(p.(pathspec{i}));
				n_def = numel(p.phi_def);
				
				if n_thisfield < n_def;
					
					if n_thisfield ~= 1
						warning(['Fewer ' pathspec{i} ' entries than phi_defs. Replicating the last entry'])
					end
					
					% Replicate last entry into needed paths
					p.(pathspec{i})( (n_thisfield+1) : n_def ) ...
						= repmat(p.(pathspec{i})(n_thisfield),[n_def-n_thisfield,1]);
					
				end
				
				% Loop over the paths, matching the segment counts
				for j = 1:n_def
					
					n_seg_def = numel(p.phi_def{j});
					n_seg_this = numel(p.(pathspec{i}){j});
					
					if n_seg_this < n_seg_def
					
						if n_seg_this ~= 1
							warning(['Fewer ' pathspec{i} ' segment entries than phi_def segments. Replicating the last entry'])
						end

						% Replicate last entry into needed paths
						p.(pathspec{i}){j}( (n_seg_this+1) : n_seg_def ) ...
							= repmat(p.(pathspec{i}){j}(n_seg_this),[n_seg_def-n_seg_this,1]);
					
					end
					
				end
	
			
				
		end
		
	end
	
end

function B = doublewrap(A,classname)
% take object A, verify that its core is of class CLASSNAME, and nest it inside
% cell functions sufficient to make it double-wrapped to return as B

	% test for one layer deep
	if iscell(A)
		
		% test for two layers deep -- make sure all the contents of A are
		% cells
		if all(cellfun(@(x) iscell(x),A)) 

			% if properly nested, do nothing
			
		else
			
			% if it's only one cell deep, add a cell layer around it
			A = {A};
			
		end
		
	else
		
		% if it's not a cell array, make it a double layered one
		A = {{A}};
		
	end

	
	% Check to make sure all contents are of the right class
	class_test = zeros(numel(A),1);
	for i = 1:numel(A)
        class_test_pieces = cellfun(@(x) isa(x,classname),A{i});
		class_test(i) = all(class_test_pieces(:));
	end

	if all(class_test)

		% If all meet the specified condition, do nothing				

	else

		% If it's not the right class, throw an error
		error(['The core of the input object should be of class ' classname ', but is of class ' class(A) '.'])

	end
	
	B = A;

	
end

