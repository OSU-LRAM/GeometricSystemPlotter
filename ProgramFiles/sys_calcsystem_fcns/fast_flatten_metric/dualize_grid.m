function [xc,yc] = dualize_grid(x,y,padding)
% Make an ndgrid whose points are at the center of the cells in a provided
% ndgrid. If padding is "true", also make a set of points that are outside
% the grid.

%%%%Get a set of x,y points at the mean position of each cell	
n_x = size(x,1);
n_y = size(x,2);

% Center-point in x direction, with robustness against non-square
% blocks (e.g., different sampling densities in different directions
xc = ((x(2:n_x,1:(n_y-1)) - x(1:(n_x-1),1:(n_y-1)))/2 + x(1:(n_x-1),1:(n_y-1)) + ...
    (x(2:n_x,2:n_y) - x(1:(n_x-1),2:n_y))/2 + x(1:(n_x-1),2:n_y))/2;

% Center-point in y direction
yc = ((y(1:(n_x-1), 2:n_y) - y(1:(n_x-1),1:(n_y-1)))/2 + y(1:(n_x-1),1:(n_y-1)) + ...
    (y(2:n_x,2:n_y) - y(2:n_x,1:(n_y-1)))/2 + y(2:n_x,1:(n_y-1)))/2;

if padding
    deltax1 = (x(2,1)-x(1,1));
    deltax2 = (x(end,1)-x(end-1,1));
    
    deltay1 = (y(1,2)-y(1,1));
    deltay2 = (y(1,end)-y(1,end-1));
    
    xc = [ xc(1,1)-deltax1 xc(1,:)-deltax1 xc(1,end)-deltax1;
           xc(:,1) xc      xc(:,end);
           xc(end,1)+deltax2 xc(end,:)+deltax2 xc(end,end)+deltax2];
       
    yc = [ yc(1,1)-deltay1 yc(1,:) yc(1,end)+deltay2;
           yc(:,1)-deltay1 yc      yc(:,end)+deltay2;
           yc(end,1)-deltay1 yc(end,:) yc(end,end)+deltay2];
       
end

end
