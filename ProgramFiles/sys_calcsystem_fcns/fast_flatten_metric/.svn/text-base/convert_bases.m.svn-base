% Find the transformation to apply to the basis functions at (x_p,y_p) in
% (x_old,y_old) to get the new basis functions at (x_p_new,y_p_new) in (x_new,y_new)
function basis_transformation = convert_bases(x_old,y_old,x_p,y_p,jacobian)

	% Interpolate the jacobian arrays
	
	basis_transformation = zeros(2,2);
	for i = 1:4
		
		basis_transformation(i) = arrayfun(@(x_p,y_p) interp2(y_old,x_old,jacobian{i},y_p,x_p,'Cubic'),x_p,y_p);
		
	end
	
	
% 	transformation = cell(2,2);
% 	for i = 1:4
% 		
% 		transformation{i} = arrayfun(@(x_p,y_p) interp2(y_old,x_old,jacobian{i},y_p,x_p,'Cubic'),x_p,y_p);
% 	
% 	end
% 	
% 	% Separate the warp arrays into metric tensors
% 	basis_transformation = cell(size(x_p));
% 
% 	for i = 1:numel(basis_transformation)
% 		
% 		% Generate a 2x2 matrix in the cell
% 		basis_transformation{i} = zeros(2,2);
% 		
% 		% Populate the matrix from the interpolated warp arrays
% 		for j = 1:4
% 			
% 			basis_transformation{i}(j) = transformation{j}(i);
% 			
% 		end
% 		
% 	end
	

end