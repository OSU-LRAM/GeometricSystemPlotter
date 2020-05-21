function jacobian = find_jacobian3(x_start,y_start,x_final,y_final,z_final)

    z_start = zeros(size(x_start));
	% Jacobian is a cell array, with each element as large as the grid
	jacobian = cell(3,3);
	
	% Populate the jacobian from old to new coordinates
	[jacobian{1,3},jacobian{1,2}, jacobian{1,1}] = gradient(x_final,z_start(1,:),y_start(1,:),x_start(:,1));
	[jacobian{2,3},jacobian{2,2}, jacobian{2,1}] = gradient(y_final,z_start(1,:),y_start(1,:),x_start(:,1));
    [jacobian{3,3},jacobian{3,2}, jacobian{3,1}] = gradient(z_final,z_start(1,:),y_start(1,:),x_start(:,1));
	

end