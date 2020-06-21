%Create grids for evaluating the connection functions
function s = create_grids(s)


    % Backwards-compatibility code for old way of specifying finite element
    % density
    if isfield(s,'finite_element_density')
        if isfield(s.density,'finite_element')
            warning('Sysf_ file has both an s.density.finite_element specification and a (deprecated) s.finite_element_density specification. The latter will be ignored, and should be removed from the sysf_ file')
        else
            warning('Sysf_ file has a (deprecated) s.finite_element_density specification. This should be replaced by an s.density.finite_element specification')
            s.density.finite_element = s.finite_element_density;
        end
    end

            
    %list of grids that will be made
    grid_list = {'vector','scalar','eval','metric_eval','metric_display','finite_element','mass_eval','coriolis_eval'};
        
    % Create a cell array to hold the grid primitives
	gridprim = cell(length(s.grid_range)/2,1);
	
    %Loop over list, creating the grids
    for i = 1:length(grid_list)
		
        % default value for any grid size
        if ~isfield(s.density,grid_list{i})
            s.density.(grid_list{i}) = 11;
        end
        
		% If necessary, duplicate out a single value given as the grid
		% density
		if numel(s.density.(grid_list{i})) == 1
			s.density.(grid_list{i}) = repmat(s.density.(grid_list{i}),length(s.grid_range)/2,1);
		end
		
		% Loop over number of dimensions
		for j = 1:length(s.grid_range)/2
			
			% Create the base vector for the grid along the jth dimension
			gridprim{j} = linspace(s.grid_range(1+(2*(j-1)))...
				,s.grid_range(2+(2*(j-1))),s.density.(grid_list{i})(j));
			
		end
        
		% Create a cell array to hold the grids
		s.grid.(grid_list{i}) = cell(length(s.grid_range)/2,1);
		
        %Grid plaid matrices
        [s.grid.(grid_list{i}){:}] = ndgrid(gridprim{:});
        
    end
    
    %%%%%%%
    % Create visualization grid if not specified elsewhere
    
    if isfield(s,'visual') && ~isfield(s.visual,'grid')
    
        % Get the number of shape variables
        n_shape = nargin(s.A_num);

        % Create a cell to hold the grid elements
        s.visual.grid = cell(n_shape,1);

        % Fill in the grid points using ndgrid. If only one set of grid points
        % is supplied, use them for all dimensions. Otherwise, use the ith set
        % of grid values for the ith shape dimension

        if isa(s.visual.grid_spacing,'numeric')
            s.visual.grid_spacing = {s.visual.grid_spacing};
        end

        [s.visual.grid{:}] = ndgrid(s.visual.grid_spacing{:});
    
    end


end