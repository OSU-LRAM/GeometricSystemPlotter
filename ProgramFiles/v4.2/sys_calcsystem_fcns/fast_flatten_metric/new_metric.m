% Produce the metric in the new coordinates, at a point specified in the old
% coordinates
function met = new_metric(x,y,M,J)

% 	% Calculate the original metric
% 	metric = arrayfun(@(x,y) M(x,y),x,y,'UniformOutput',false);
	
	% Calculate the inverse jacobian of the transformation to the grid
	jacobian = (J(x,y));%arrayfun(@(x,y) W(x,y),x,y,'UniformOutput',false);
	inv_jacobian = inv(jacobian);
	
	% Multiply the convert the metric
% 	met = cellfun(@(j,m) j' * m * j,metric,jacobian,'UniformOutput',false);
	met = inv_jacobian' * M(x,y) * inv_jacobian;
end