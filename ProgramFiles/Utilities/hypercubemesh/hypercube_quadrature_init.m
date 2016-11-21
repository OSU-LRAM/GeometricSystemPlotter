function [locations, weights] = hypercube_quadrature_init(N)
% Generate the quadrature point locations and weights for an N-dimensional
% hypercubic finite element

	% For now, assume third-order quadrature
	
	%%%%%
	% Get the locations of the points (as per Becker,Carey,Oden)
	
	% Cell to hold the quadrature point locations
	locations = cell(N,1);
	
	% coordinates of the quadrature points along each axis, replicated as a
	% cell along each dimension
	point_coordinates = repmat({sqrt(.6)*[-1 0 1]'},[N 1]);
	
	% Use ndgrid to generate the mesh of quadrature locations
	if N~=1
		[locations{:}] = ndgrid(point_coordinates{:});
	else
		locations = point_coordinates;
	end
	
	% flatten the location arrays into lists
	locations = cellfun(@(x) x(:),locations,'UniformOutput',false);
	
	%%%%%%%%
	% Get the weights
	
	% cell array to hold the weight contributions
	weight_contribution = cell(N,1);
	
	% weight basis along each axis, replicated as a cell array
	weight_basis = repmat({[5 8 5]'/9},[N 1]);
	
	% use ndgrid to generate the mesh of quadrature weights
	if N~=1
		[weight_contribution{:}] = ndgrid(weight_basis{:});
	else
		weight_contribution = weight_basis;
	end
	
	% Concatenate the weight contributions along the N+1 dimension, then
	% multiply down to get each weight
	weights = prod( cat(N+1,weight_contribution{:}), N+1);
	
	% flatten the weights into a list
	weights = weights(:);
	
	
	
end