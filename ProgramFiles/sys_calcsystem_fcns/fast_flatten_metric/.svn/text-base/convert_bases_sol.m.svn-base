% Find the transformation to apply to the basis functions at (x_p,y_p) in
% (x_old,y_old) to get the new basis functions at (x_p_new,y_p_new) in (x_new,y_new)
function basis_transformation = convert_bases_sol(x_old,y_old,x_p,y_p,sol,t)

	points = deval(sol,t);
	
	[x_t, y_t] = deal(zeros(size(x_old)));
	x_t(:) = points(1:(end/2));
	y_t(:) = points((1+end/2):end);
	
	% Get the new Jacobian
	jacobian = find_jacobian(x_old,y_old,x_t,y_t);
	
	% Interpolate the jacobian arrays
	
	basis_transformation = zeros(2,2);
	for i = 1:4
		
		basis_transformation(i) = arrayfun(@(x_p,y_p) interp2(y_old,x_old,jacobian{i},y_p,x_p,'Cubic'),x_p,y_p);
		
	end
	
end