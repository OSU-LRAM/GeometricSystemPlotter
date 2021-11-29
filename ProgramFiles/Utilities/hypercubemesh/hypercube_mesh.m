% generate a hypercubic mesh over an arbitrarily-dimensional space

% inputs:
    % grid: ndgrid
% outputs:
    % nodes:
        % list (matrix) of points. each row contains columns of
        % [x, y, z, w, ...]
    % cubes:
        % list (matrix) of points in each cube. each row contains columns
        % of [corner1, corner2, ...] points

function [nodes,cubes] = hypercube_mesh(grid)

    %% generate node information
	% turn grid into node list
	nodes = grid_to_columns(grid);
    
	% node list metadata
	node_nums = zeros(size(grid{1}));
	node_nums(:) = 1:numel(node_nums); %indexing of nodes
	nodes_per_cube = 2^numel(grid); %2^N nodes per cube
	l = size(grid{1}); %lengths of dimensions
	
	%% build list of nodes in each hypercube
	
	% prime cubes array (one fewer cube than node in each direction)
	cubes = zeros(prod(size(node_nums)-1),nodes_per_cube);
	n_cubes = size(cubes,1);
    
    % get possible ways to vary start/end points (done with binary)
    bin_strs = dec2bin(0:nodes_per_cube - 1);
    bin_mat = zeros(size(bin_strs));
    bin_mat(:) = arrayfun(@str2num, bin_strs(:));
    % use binary values to assign cube vertices
    for vertex = 1:nodes_per_cube
        % get nodes (in each dimension)
        idxs = cell(size(l));
        for dim = 1:length(l)
            offset = bin_mat(vertex, dim);
            idxs{dim} = 1 + offset : l(dim) - 1 + offset;
        end
        % get the id of each node stored by matlab (rather than coordinate)
        cubes(:, vertex) = reshape(node_nums(idxs{:}), n_cubes, 1);
    end
end