function [h, J] = backbone_from_general_curvature(curvdef,cparams,L)

% Set the integration limits. Internally, everything uses [-0.5 0.5], and
% then scales by the length factor at the end.
all_limits = [-0.5 0 0.5];

% Get backbone theta as a function of s
theta_fun = curvdef(cparams,'angle');

% Get the x and y locations of points on the backbone
locus_sol = ode_multistart(@ode45, @(s,h) locus_helper(s,theta_fun),all_limits,0,[0;0]);

% Return the locus data as 3xn matrix with x,y,theta at each of n points
h_norm = @(s) [locus_sol(torow(s)); theta_fun(torow(s))]; % in normalized internal coordinates
h = @(s) [L*locus_sol(torow(s)/L); theta_fun(torow(s)/L)]; % with actual length

%%%%%%%%%%%%%%
% Get the jacobian to body point velocities

dcurvdef = curvdef(cparams,'dcurvature_int');

%Jacobian of theta function is the derivative of the curvature with respect
%to each of the parameters, as definded in the curvdef function
J_theta_fun = dcurvdef; %ode_multistart(@ode45,@(s,J) dcurvdef(s),all_limits,0,zeros(length(dcurvdef(0)),1));

% SE(2) adjoint integration to get the x and y jacobians
jacobian_sol = ode_multistart(@ode45,@(s,J) J_helper(s,J,dcurvdef,h_norm(s)),all_limits,0,zeros(2*length(dcurvdef(0)),1));

% Concatenate xy and theta jacobians.
J = @(s) cat(1,reshape(L*jacobian_sol(toz(s/L)),2,[],length(s)),ipermute(J_theta_fun(s/L),[2,3,1]));



end


% Get the derivative of the locus along the backbone length
function dh = locus_helper(s,theta_fun)

	% Build the rotation matrix
	theta = theta_fun(s);
	R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
		
	% Rotate the basic tangent vector to point along the backbone
	dh = R * [1; 0];
	
end

% get the derivative of the jacobian along the backbone length (note that J
% is only here for ODE45, and isn't actually used)
function dJ = J_helper(s,J,dcurvdef,h) %#ok<INUSL>

    nvars = length(dcurvdef(0));

	dJ = zeros(2,nvars);
    
    dcurvdef_eval = dcurvdef(s);
	
	for i = 1:nvars
		
		int_basis = dcurvdef_eval(i,:);
		
		dJ(:,i) = [-sin(h(3)) -cos(h(3));
			cos(h(3)) -sin(h(3))] * [int_basis; 0];
		
	end
	
	%Reshape into a vector for ODE
	dJ = dJ(:);

end

function X = torow( X ) % Make the argument a row

	X = X(:).';
    
end

function Y = toz( X )

	%first, make X a row
	X = torow(X);
	
	% then, make it a 3rd-dimension column
	Y(1,1,:) = X;
    
end