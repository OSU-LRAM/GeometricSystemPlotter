function [V,grid] = plotting_interp(V,grid,resolution,resolution_choice)

	gridres = resolution.(resolution_choice);
	gridrange = resolution.([resolution_choice '_range']);


	% If necessary, duplicate out a single value given as the grid
	% density
	if numel(gridres) == 1
		gridres = repmat(gridres,numel(grid),1);
	end

	% Loop over number of dimensions
    gridprim = cell(numel(grid),1);
	for i = 1:numel(grid)

		% Create the base vector for the grid along the jth dimension
		gridprim{i} = linspace(gridrange(1+(2*(i-1)))...
			,gridrange(2+(2*(i-1))),gridres(i));

	end

	% Create a cell array to hold the grids
	plot_grid = cell(numel(grid),1);

	%Grid plaid matrices
	[plot_grid{:}] = ndgrid(gridprim{:});
	
	%%
	% Interpolate vector fields onto plotting grid
	for i = 1:numel(V)
		V{i} = interpn(grid{:},V{i},plot_grid{:},'cubic');
	end
	
	%%
	% Replace the grid with the plot grid
	grid = plot_grid;
	
end