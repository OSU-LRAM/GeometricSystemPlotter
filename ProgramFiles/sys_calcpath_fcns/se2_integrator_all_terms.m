% Function to evaluate displacement and differential cost at each time
function V = se2_integrator_all_terms(t,X,s,phi_fun,dphi_fun)

	% State order is bvi_orig,bvi_opt,disp_orig,disp_opt

	% Get the shape and shape derivative at the current time
	shape = phi_fun(t);
    shape = shape(:);
	shapelist = num2cell(shape);
	dshape = dphi_fun(t);
    dshape=dshape(:);
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
    
	% Get the local connection at the current time, in both sets of
	% coordinates
	A_original = cellfun(@(Y) -interpn(s.grid.eval{:},Y,shapelist{:},'spline'),s.vecfield.eval.content.Avec);

	A_optimized = cellfun(@(Y) -interpn(s.grid.eval{:},Y,shapelist{:},'spline'),s.vecfield.eval.content.Avec_optimized);
			
	
	% Get the body velocity at the current time
	xi_original = - A_original * dshape(:);
	xi_optimized = - A_optimized * dshape(:);
	
	% Rotate body velocity into world frame for displacement integrals
	theta_original = X(9);
	TeLg_original = [cos(theta_original) -sin(theta_original) 0;...
		sin(theta_original) cos(theta_original) 0; 0 0 1];
	
	theta_optimized = X(12);
	TeLg_optimized = [cos(theta_optimized) -sin(theta_optimized) 0;...
		sin(theta_optimized) cos(theta_optimized) 0; 0 0 1];

	% Output the velocities
	V = [xi_original;xi_optimized;TeLg_original*xi_original;TeLg_optimized*xi_optimized];
	


end