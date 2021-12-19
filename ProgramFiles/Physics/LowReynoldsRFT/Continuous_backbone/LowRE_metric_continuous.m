function Mp = LowRE_metric_continuous(geometry,physics,shapeparams)
% Calculate the dissipation power metric for a set of curvature bases

	% Specified integration limits
	int_limit = geometry.length*[-0.5 0.5];
	
	% Define the tangential, lateral drag matrix for unit/double drag
	drag_matrix = [1 0; 0 physics.drag_ratio]*physics.drag_coefficient;

	% Get the backbone locus, Jacobian, and Local Connection functions
    [A,h,J] = LowRE_local_connection(geometry,physics,shapeparams);

	% Integrate along the body for the power metric
	Mp_sol = ode45(@(s,Mp) dMetric(s,Mp,A,h(s/geometry.length),J(s/geometry.length),drag_matrix),int_limit,zeros(length(shapeparams)));%zeros(length(shapeparams)^2,1));

    % Form the components of the metric into the standard NxN symmetric matrix form
	Mp = reshape(deval(Mp_sol,int_limit(end)),length(shapeparams),[]);

end

function dMp = dMetric(s,Mp,A,h,J,drag_matrix) %#ok<INUSL>
% Calculate the local contribution to the power metric

    % Calculate the jacobian from shape variables to body velocity of this
    % point on the system:
	
    % Body velocity for a point on the robot is sum of right-translated
    % system body velocity and velocity within body frame, left translated
    % into local coordinates
    A_point = -(TgLginv(h) * (-TeRg(h)*A + J));
    
    % the xy jacobian from shape variables to body velocity of points on the system is the
    % first two rows of the local connection.
    localJ = -A_point(1:2,:);

    % Contribution to the metric from this point is its local drag matrix,
    % pulled back into joint coordinates by its Jacobian.
	dMp = localJ.'*drag_matrix*localJ;
	
    % Convert metric contribution into column vector for ODE integration
	dMp = dMp(:);
	
end