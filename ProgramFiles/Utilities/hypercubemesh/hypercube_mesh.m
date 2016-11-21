function [nodes,cubes] = hypercube_mesh(grid)
% Generate a normalized tetrahedral mesh over a hypercubic region with a
% given density

	% turn grid into node list
	nodes = grid_to_columns(grid);

	% Number the nodes
	node_nums = zeros(size(grid{1}));
	node_nums(:) = 1:numel(node_nums);
	
	% 2^N nodes per cube
	nodes_per_cube = 2^numel(grid);
	
	% number of nodes in each direction
	l = size(grid{1});
	
	% Build the list of nodes in each hypercube
	
	% Prime cubes array (one fewer cube than node in each direction
	cubes = zeros(prod(size(node_nums)-1),nodes_per_cube);
	
	% number of cubes
	n_cubes = size(cubes,1);
			
	% Expand algorithm for higher dimensions
	switch (numel(grid))
		
		case 2
			
			cubes(:,1) = reshape(node_nums(1:l(1)-1,1:l(2)-1),n_cubes,1);
			cubes(:,2) = reshape(node_nums(2:l(1),1:l(2)-1),n_cubes,1);
			cubes(:,3) = reshape(node_nums(2:l(1),2:l(2)),n_cubes,1);
			cubes(:,4) = reshape(node_nums(1:l(1)-1,2:l(2)),n_cubes,1);
			
		case 3
			
			cubes(:,1) = reshape(node_nums(1:l(1)-1,1:l(2)-1,1:l(3)-1),n_cubes,1);
			cubes(:,2) = reshape(node_nums(2:l(1),1:l(2)-1,1:l(3)-1),n_cubes,1);
			cubes(:,3) = reshape(node_nums(2:l(1),2:l(2),1:l(3)-1),n_cubes,1);
			cubes(:,4) = reshape(node_nums(1:l(1)-1,2:l(2),1:l(3)-1),n_cubes,1);
			
			cubes(:,5) = reshape(node_nums(1:l(1)-1,1:l(2)-1,2:l(3)),n_cubes,1);
			cubes(:,6) = reshape(node_nums(2:l(1),1:l(2)-1,2:l(3)),n_cubes,1);
			cubes(:,7) = reshape(node_nums(2:l(1),2:l(2),2:l(3)),n_cubes,1);
			cubes(:,8) = reshape(node_nums(1:l(1)-1,2:l(2),2:l(3)),n_cubes,1);
			
		otherwise
			
			error('hypercube_mesh not yet implemented for N>3 dimensions');
			
	end			
			

end