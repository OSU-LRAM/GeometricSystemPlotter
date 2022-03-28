% Function to evaluate velocity and differential cost at each time for ODE
% solver
function dX = helper_function(t,X,s,gait,ConnectionEval)

	% X is the accrued displacement and cost

	[xi, dcost] = get_velocities(t,s,gait,ConnectionEval);
		
	% Rotate body velocity into world frame
	theta = X(3);
	v = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*xi;
		
	% Combine the output
	dX = [v;dcost];
	

end