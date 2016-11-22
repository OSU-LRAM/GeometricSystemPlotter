% Turn a cell array of grids into a simple array of columns
function columns = grid_to_columns(grid)

	columns = cell2mat(arrayfun(@(g) g{1}(:),grid(:)','UniformOutput',false));
	
end