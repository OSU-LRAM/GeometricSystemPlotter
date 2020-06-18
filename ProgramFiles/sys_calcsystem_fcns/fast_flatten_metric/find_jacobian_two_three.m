function jacobian = find_jacobian_two_three(x_start,y_start,x_final,y_final,z_final)

	% Jacobian is a cell array, with each element as large as the grid
	jacobian = cell(3,2);
	
	% Populate the jacobian from old to new coordinates
	[jacobian{1,2}, jacobian{1,1}] = gradient(x_final,y_start(1,:),x_start(:,1));
	[jacobian{2,2}, jacobian{2,1}] = gradient(y_final,y_start(1,:),x_start(:,1));
	[jacobian{3,2}, jacobian{3,1}] = gradient(z_final,y_start(1,:),x_start(:,1));	

end