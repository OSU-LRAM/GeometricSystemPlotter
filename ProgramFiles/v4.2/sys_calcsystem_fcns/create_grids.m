%Create grids for evaluating the connection functions
function s = create_grids(s)

    %list of grids that will be made
    grid_list = {'vector','scalar','eval','metric_eval'};
        
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

end