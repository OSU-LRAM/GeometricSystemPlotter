function gridSpacing = extractGridSpacing(grid)
% Extract grid spacing vectors from an ndgrid

% create a cell with as many entries as there are dimensions in the space
gridSpacing = cell(size(grid));

% Make a cell array with a colon in the first entry, and padded out with
% ones to the dimension of the space
cellones = num2cell(ones(numel(grid)-1,1));
extractionvector = [':',torow(cellones)];

% Loop over the grid
for idx = 1:numel(grid)
    
    % Make a list of all indices other than the current dimension
    otherindices = 1:numel(grid);
    otherindices(otherindices==idx) = [];
    
    % Permute the grid so that the current index is the first dimension
    permutedgrid = permute(grid{idx},[idx,otherindices]);
    
    % Pull out the first column of the permuted grid
    gridSpacing{idx} = permutedgrid(extractionvector{:});
end

end