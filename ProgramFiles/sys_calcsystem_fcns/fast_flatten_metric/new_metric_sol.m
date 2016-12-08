% Produce the metric in the new coordinates, at a point specified in the old
% coordinates
function met = new_metric_sol(x,y,t,M,J_fun)

	% Calculate the inverse jacobian of the transformation to the grid
	jacobian = (J_fun(x,y,t));%arrayfun(@(x,y) W(x,y),x,y,'UniformOutput',false);
	inv_jacobian = inv(jacobian);
	
	% Multiply the convert the metric
% 	met = cellfun(@(j,m) j' * m * j,metric,jacobian,'UniformOutput',false);
	met = inv_jacobian' * M(x,y) * inv_jacobian;
end