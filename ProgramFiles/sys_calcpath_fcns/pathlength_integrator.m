function dS = pathlength_integrator(t,y,s,phi_fun,dphi_fun) %#ok<INUSL>

	% Get the shape and shape derivative at the current time
	shape = phi_fun(t);
	shapelist = num2cell(shape);
	dshape = dphi_fun(t);
	
	% Evaluate the metric tensor

	M = cellfun(@(Y) interpn(s.grid.metric_eval{:},Y,shapelist{:}),s.metricfield.metric_eval.content.metric);
	
	% get the contribution to the pathlength
	dS = sqrt(dshape(:)' * M * dshape(:));
	

	
end