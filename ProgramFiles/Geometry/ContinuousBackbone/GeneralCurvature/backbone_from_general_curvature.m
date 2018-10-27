function [h, J, J_full] = backbone_from_general_curvature(curvdef,cparams,L)
% Integrate backbone geometry from a generalized curvature definition
% function and a set of parameters on this curvature function. The backbone
% is of total length L, placed such that its midpoint is at the origin and
% its midpoint tangent is aligned with the x axis.
% 
%
% Inputs: 
%
% curvdef: function handle to a curv_ function whose m-file is within the
%   current matlab path (for sysplotter, this will most likely be in the
%   user's Systems/ directory, near the sysf file for the system. See
%   make_curvdef for how to create this mfile.
%
% cparams: Vector of parameters to be passed to the curvdef file
%
% L: The actual length of the backbone, in "real world" units. This value
%   scales the geometry of the integrated backbone shape
%
% Outputs:
%
% h is a function on the domain s = [-.5 .5] (normalized unit-length backbone) that returns the position and
%   orientation (in world-units) of the frame on the backbone that is at point s along the
%   backbone (in backbone units)
% J is the Jacobian of h, so that its product with a shape velocity gives
%   the rate of change of h -- the "world velocity" of the backbone point
%   taking the midpoint-tangent frame as "fixed"
% J_full is the Jacobian mapping from system-body-velocity and shape
%   velocity to the body velocity (relative to a non-moving frame) of the
%   point on the body


% Set the integration limits. Internally, everything uses [-0.5 0.5], and
% then scales by the length factor at the end.
all_limits = [-0.5 0 0.5];

% Get backbone theta as a function of s
theta_fun = curvdef(cparams,'orientation');

% Get the x and y locations of points on the backbone
locus_sol = ode_multistart(@ode45, @(s,h) locus_helper(s,theta_fun),all_limits,0,[0;0]);

% Return the locus data as 3xn matrix with x,y,theta at each of n points
h_norm = @(s) [locus_sol(torow(s)); theta_fun(torow(s))]; % in normalized internal coordinates
h = @(s) [L*locus_sol(torow(s)); theta_fun(torow(s))]; % with actual length

%%%%%%%%%%%%%%
% Get the jacobian to body point velocities

if nargout > 1

    dcurvdef = curvdef(cparams,'dcurvature_int');

    %Jacobian of theta function is the derivative of the curvature with respect
    %to each of the parameters, as definded in the curvdef function
    J_theta_fun = dcurvdef; %ode_multistart(@ode45,@(s,J) dcurvdef(s),all_limits,0,zeros(length(dcurvdef(0)),1));

    % SE(2) adjoint integration to get the x and y jacobians
    jacobian_sol = ode_multistart(@ode45,@(s,J) J_helper(s,J,dcurvdef,h_norm(s)),all_limits,0,zeros(2*length(dcurvdef(0)),1));

    % Concatenate xy and theta jacobians.
    J = @(s) cat(1,reshape(L*jacobian_sol(toz(s)),2,[],length(s)),permute(J_theta_fun(s),[3,2,1]));

end

%%%%%%%
% Get the full Jacobian from system body velocity and shape velocity to
% body velocity of the point on the body (relative to fixed world frame)
if nargout > 2
    
    J_full = @(s) [Adjinv(h(s)) TgLginv(h(s))*J(s)];
    
end


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
		
		int_basis = dcurvdef_eval(:,i);
		
		dJ(:,i) = [-sin(h(3)) -cos(h(3));
			cos(h(3)) -sin(h(3))] * [int_basis; 0];
		
	end
	
	%Reshape into a vector for ODE
	dJ = dJ(:);

end



