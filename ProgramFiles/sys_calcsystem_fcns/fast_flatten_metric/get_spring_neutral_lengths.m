% This finds the neutral length of each spring, which is its length if the
% box it's part of is premultiplied by sqrt(M)^-1. For springs that are part of
% two boxes, the neutral length is taken as the average across the two
% boxes
function [metric_lengths,mean_length] = get_spring_neutral_lengths(springs,blocks,start_deltas,M)

% 	%%%%Get a set of x,y points at the mean position of each cell	
% 	n_x = size(x,1);
% 	n_y = size(x,2);
% 	
% 	% Center-point in x direction, with robustness against non-square
% 	% blocks (e.g., different sampling densities in different directions
% 	xc = ((x(2:n_x,1:(n_y-1)) - x(1:(n_x-1),1:(n_y-1)))/2 + x(1:(n_x-1),1:(n_y-1)) + ...
% 		(x(2:n_x,2:n_y) - x(1:(n_x-1),2:n_y))/2 + x(1:(n_x-1),2:n_y))/2;
% 	
% 	% Center-point in y direction
% 	yc = ((y(1:(n_x-1), 2:n_y) - y(1:(n_x-1),1:(n_y-1)))/2 + y(1:(n_x-1),1:(n_y-1)) + ...
% 		(y(2:n_x,2:n_y) - y(2:n_x,1:(n_y-1)))/2 + y(2:n_x,1:(n_y-1)))/2;
% 	
% 	% Evaluate the metric at the center of each cell,
% 	metric = celltensorconvert(...
% 		cellfun(@(m) interpn(x,y,m,xc,yc,'cubic'),M,'UniformOutput',false));	

    metric = celltensorconvert(M);
	
	%%%%%%%%%%
	% Make a cell array of the coordinates of the four nodes associated
	% with each block, translated to be centered on (xc,yc), and then
	% transformed by the inverse of the metric
	
	% make a cell array with block data in it
	blocks_cell = num2cell(blocks,2);
    
    % evaluate the length of each spring in the block by multiplying the
    % spring deltas in the block by the metric, then taking the square root
    % of the diagonal elements of the resulting array (could speed this up
    % a bit by evaluating only the diagonal elements, this is more concise)
    spring_lengths_cell = cellfun(@(block,metric) sqrt(diag(start_deltas(block,:)*metric*start_deltas(block,:)'))'...
        ,blocks_cell,metric(:),'UniformOutput',false)';
    
    % Flatten the spring lengths into a single matrix
    spring_lengths = vertcat(spring_lengths_cell{:});
    
    % Average the lengths of springs that are in multiple cells, using
    % arithmetic mean
	metric_lengths = zeros(size(springs,1),1);
	for i = 1:length(metric_lengths)
		
		metric_lengths(i) = mean(spring_lengths(blocks == i));
		
    end
    
    % Calculate the geometric mean of the spring lengths
    mean_length = geomean(metric_lengths);
    
end