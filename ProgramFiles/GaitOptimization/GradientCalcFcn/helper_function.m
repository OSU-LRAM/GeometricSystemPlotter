% Function to evaluate velocity and differential cost at each time for ODE
% solver
function dX = helper_function(t,X,s,gait,ConnectionEval)

    for dim = 1:numel(gait.phi_def)
        samplePoints{dim} = gait.phi_def{dim}(t);
    end

    As = cellfun(@(C) -interpn(s.grid.eval{:},C,...
        samplePoints{:},'spline'),s.vecfield.eval.content.Avec,...
        'UniformOutput',false);
    As = celltensorconvert(As);
    As = As{:};

    %Otherwise, we're not doing inertial so don't need dM, but we
    %do care about the metric and connection
    metrics = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
        samplePoints{:},'spline'),s.metricfield.metric_eval.content.metric,...
        'UniformOutput',false);
    metrics = celltensorconvert(metrics);
    metrics = metrics{:};
    
	% X is the accrued displacement and cost

	[xi, dcost] = get_velocities(t,s,gait,ConnectionEval,...
                As,metrics,1);
	% Rotate body velocity into world frame
	theta = X(3);
	v = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*xi;
		
	% Combine the output
	dX = [v;dcost];
	

end