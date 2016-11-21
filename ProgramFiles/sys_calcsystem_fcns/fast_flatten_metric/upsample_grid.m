function [x_up, y_up] = upsample_grid(x,y,up_factor)
% Increase the density of a grid by UP_FACTOR. This code uses the "natural"
% interpretation of increasing the density, where an UP_FACTOR of 2 places
% a new gridline in between each existing gridline, rather than simply
% multiplying the existing number of gridlines by 2. New grid lines are
% linearly-interpolated between existing grid lines.

% Extract the generating vectors for the grid
x_gen = x(:,1)';
y_gen = y(1,:);

% % Get the lengths of the generating vectors
% x_size = length(x_gen);
% y_size = length(y_gen);

% Get the differences along the vectors
x_diff = x_gen(2:end)-x_gen(1:end-1);
y_diff = y_gen(2:end)-y_gen(1:end-1);

% Stack the vectors with interpolated values (first value in each vector is
% skipped here, but restored later)
x_stacker = zeros(numel(x_diff),up_factor);
y_stacker = zeros(numel(y_diff),up_factor);

for i = 1:up_factor
	
	x_stacker(i,:) = x_gen(1:end-1) + i/up_factor*x_diff;
	y_stacker(i,:) = y_gen(1:end-1) + i/up_factor*y_diff;
	
end

% Turn the interpolated numbers into single vector, and restore the first
% value in each
x_up_gen = [x_gen(1);x_stacker(:)];
y_up_gen = [y_gen(1);y_stacker(:)];

% Re-plaid the vectors
[x_up, y_up] = ndgrid(x_up_gen,y_up_gen);



% % Start by getting the dimensions of the initial grid
% gridsize = size(x);
% 
% % Make a pair of new arrays to hold the expanded grid
% newgridsize = (up_factor*(gridsize-1)+1);
% x_up = zeros(newgridsize);
% y_up = x_up;
% 
% 
% % Place interpolated values into the new grids
% diffx = x(2:end,:)-x(1:end-1,:);
% diffy = y(:,2:end)-y(:,1:end-1);
% for i = 0:up_factor;
% 	
% 	x_up(i+(1:up_factor:(end-1)),:) = x(1:(end-1),:) + (i/*diffx;
% 	y_up(:,i+(1:up_factor:(end-1))) = y(:,1:(end-1)) + i*diffy;
% 	
% end
