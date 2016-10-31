function corners = find_corners(minvals,maxvals)
% Find the corners of a hypercube extending from minvals to maxvals

	% Get the number of dimensions
	n_dim = length(minvals);
	
	
	% Preallocate the corners
	corners = zeros(2^(n_dim),n_dim);
	
	% Set up selection vector
	selection = zeros(1,n_dim);
	
	% Build the corner vectors
	for i = 1:2^(n_dim)
		
		
		corners(i,:) = minvals.*selection+maxvals.*(~selection);
		
		% Cycle the selection vector
		for j = 1:n_dim
			
			if mod(i,2^(j-1)) == 0;
				
				selection(j) = ~selection(j);
				
			end
			
		end
		
	end

end