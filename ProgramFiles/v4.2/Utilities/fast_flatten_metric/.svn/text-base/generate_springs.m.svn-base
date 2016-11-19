function [springs, blocks] = generate_springs(x)

	%number of nodes in total, x, and y directions
	n_t = numel(x);
	n_x = size(x,1);
	n_y = size(x,2);

	%nodes to which the vertical springs are attached
	vertical_springs = [(1:(n_t-n_x))' (1:(n_t-n_x))'+n_x];
	
	%nodes to which the horizontal springs are attached
	all_nodes = (1:n_t)';
	rightmost = n_x:n_x:n_t;
	not_rightmost = all_nodes;
	not_rightmost(rightmost) = [];
	horizontal_springs = [not_rightmost,not_rightmost+1];
	
	%nodes to which the +45 degree springs are attached
	topmost = (n_t-n_x+1):n_t;
	not_right_nor_topmost = all_nodes;
	not_right_nor_topmost([rightmost topmost]) = [];
	p45_springs = [not_right_nor_topmost, not_right_nor_topmost+1+n_x];
	
	%nodes to which the -45 degree springs are attached
	bottommost = 1:n_x;
	not_right_nor_bottommost = all_nodes;
	not_right_nor_bottommost([rightmost bottommost]) = [];
	m45_springs = [not_right_nor_bottommost, not_right_nor_bottommost+1-n_x];
	
	% Group the spring nodes, tagging the different blocks so as not to
	% lose the data
	springs = [vertical_springs ones(size(vertical_springs,1),1);
		horizontal_springs 2*ones(size(horizontal_springs,1),1);
		p45_springs 3*ones(size(p45_springs,1),1);
		m45_springs 4*ones(size(m45_springs,1),1)];
	
	
	%%%%%%%%%
	% blocks
	
	% Initialize an array with the number of cells between the nodes, and
	% space for the six springs in that cell
	blocks = zeros((n_x-1)*(n_y-1),6);

	%Numeration of the blocks
	blocklist = (1:size(blocks,1))';
	
	
	% The left edge of the block
	blocks(:,1) = blocklist+floor((blocklist-1)/(n_x-1));
	
	% The right edge of the block is the vertical spring with number offset
	% by 1 from the left spring
	blocks(:,2) = blocks(:,1)+1;
	
	% The bottom edge of the block is the horizontal spring with the same 
	% number as the block. Horizontal spring numbers are offset by the
	% number of vertical springs
	blocks(:,3) = blocklist+size(vertical_springs,1);

	% The top edge of the block is the horizontal spring numer offset
	% by n_x-1. Horizontal spring numbers are offset by the
	% number of vertical springs
	blocks(:,4) = blocklist+(n_x-1)+size(vertical_springs,1);
	
	% The diagonal springs are numbered the same as the blocks. p45 springs
	% are offset by the number of vertical and horizontal springs
	blocks(:,5) = blocklist+size(vertical_springs,1)+size(horizontal_springs,1);
	blocks(:,6) = blocklist+size(vertical_springs,1)+size(horizontal_springs,1)+size(p45_springs,1);


end