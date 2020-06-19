% Take points in the one set of coordinates, and express them in the
% power_normalized coordinates
function [x_new, y_new, z_new] = convert_points(x_start,y_start,x_end,y_end,x_p,y_p)

	%Identify any NaN elements
	nanIDx = isnan(x_p);
	nanIDy = isnan(y_p);
	
	%Replace the NaN with valid holding values
	x_p(nanIDx) = min(x_start(:));
	y_p(nanIDy) = min(y_start(:));
	
	% odd ordering of arguments is ndgrid/meshgrid incompatibility
	x_new = interpn(x_start,y_start,x_end,x_p,y_p,'cubic');
	y_new = interpn(x_start,y_start,y_end,x_p,y_p,'cubic');

	%Reinsert the NaNs
	x_new(nanIDx) = NaN;
	y_new(nanIDy) = NaN;
    
    
    z_new = zeros(size(x_new));
	
	
end