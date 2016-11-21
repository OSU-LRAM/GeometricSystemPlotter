function [x_new, y_new] = convert_points_sol_interp(x_start,y_start,sol,t,x_old,y_old)
% Convert points at x_old,y_old to their locations at time t in the
% relaxation process

	points = deval(sol,t);
	
	[x_t, y_t] = deal(zeros(size(x_start)));
	x_t(:) = points(1:(end/2));
	y_t(:) = points((1+end/2):end);
	
	[x_new, y_new] = convert_points(x_start,y_start,x_t,y_t,x_old,y_old);

end

