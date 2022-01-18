function dS = pathlength_integrator(t,y,s,phi_fun,dphi_fun) %#ok<INUSL>

	% Get the shape and shape derivative at the current time
	shape = phi_fun(t);
    shape = shape(:);
	shapelist = num2cell(shape);
	dshape = dphi_fun(t);
    dshape = dshape(:);
    n_dim=length(s.vecfield.eval.content.Avec_optimized(1,:));
    
    if length(shape)<n_dim
        shape=[shape;zeros(n_dim-length(shape),1)];
    end
    shapelist = num2cell(shape);
    if length(dshape)<n_dim
        dshape=[dshape;zeros(n_dim-length(dshape),1)];
    end
    
    if length(shape)>n_dim
        shape=shape(1:n_dim);
    end
    shapelist = num2cell(shape);
    if length(dshape)>n_dim
        dshape=dshape(1:n_dim);
    end    
    
	% Evaluate the metric tensor

	M = cellfun(@(Y) interpn(s.grid.metric_eval{:},Y,shapelist{:}),s.metricfield.metric_eval.content.metric);
	
	% get the contribution to the pathlength
	dS = sqrt(dshape(:)' * M * dshape(:));
	

	
end