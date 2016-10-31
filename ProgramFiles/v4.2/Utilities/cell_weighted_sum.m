function output_sum = cell_weighted_sum(weights,cells)
% Build a new function that is the weighted sum of the functions in the
% cell array funs, with weights given by weights

	% Multiply each cell by its weight
	wcells = cellfun(@(x,y) x*y, num2cell(weights(:)),cells(:),'UniformOutput',false);
	
	% Concatenate the cells across a one-higher dimension
	n_dim = size(weights,1);
	catwcells = cat(n_dim+1,wcells{:});
	
	if n_dim ~= 1
		output_sum = sum(catwcells,n_dim+1);
	else
		output_sum = catwcells;
	end
	
end