function weight = weight_away_from_singularities(singularity_locs, grid)

		% Get the number of dimensions in the vector fields
		n_dim = numel(size(grid{1}));

		% Find all singularities
		singular_indices = max(cat(n_dim+1,singularity_locs{:}),[],n_dim+1);
		
		% Get the location of all singularities
		% Turn a cell array of grids into a simple array of columns
		all_points = grid_to_columns(grid);
		singular_points = all_points(singular_indices(:)==1,:);
				
		% Get the distance from each point to the each singularity
		singular_distance = pdist2(all_points,singular_points);
		
		% Get the shortest distance from each point to a singularity
		min_singular_distance = min(singular_distance,[],2);
		
		% Reshape this back into the original grid
		min_singular_distance_gridded = reshape(min_singular_distance,size(grid{1}));
		
		% Set the weight to the distance from the singularity
		%weight = (min_singular_distance_gridded).^2;
         weight = abs(atan(min_singular_distance_gridded));
		
		% Add in a fudge factor to avoid zero-values at singularities
		weight = weight+0.01*max(weight(:));
        
       
       
		
end