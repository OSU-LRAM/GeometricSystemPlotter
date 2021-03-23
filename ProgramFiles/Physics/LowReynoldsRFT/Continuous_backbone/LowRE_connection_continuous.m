function [A, h, J,J_full,Omega]= LowRE_connection_continuous(geometry,physics,shapeparams)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
% geometry: structure defining system geometry
%      geometry.type: how the system geometry is defined 
%         (e.g., links, curvature basis, general curvature)
%      geometry.function: map from shape variables to local backbone
%         deformation (e.g., curvature or joint angles)
%      geometry.length: total length of swimmer
% cparams: value of shape variables

% c: drag per unit length
% drag_ratio: ratio of lateral to longitudinal drag



%Generate backbone geometry and Jacobian from its local definition
[h,J,J_full] = backbone(geometry,shapeparams);

% Itegrate from one halflength before the midpoint to one halflength after it
int_limit = [-0.5 0.5];

% Now integrate to get the pfaffian
Omega_sol = ode45( @(s,Omega) LowRE_Pfaffian_infinitesimal(s,h(s),J(s),J_full(s),geometry.length,physics.drag_coefficient,physics.drag_ratio),int_limit,zeros(3,3+length(shapeparams)));

% Reshape the terms of the Pfaffian into a matrix of the correct dimension
Omega = reshape(deval(Omega_sol,int_limit(end)),3,[]);

% Calculate the local connection by multiplying the inverse of the first
% block in the Pfaffian by the second block
A = Omega(:,1:3)\Omega(:,4:end);


end


function dOmega = LowRE_Pfaffian_infinitesimal(s,h,J,J_full,lambda,c,drag_ratio) %#ok<INUSL>
% Calculate the derivative of the local connection as it's built up along
% the backbone

	% Convert velocity to local velocity
	gdot_to_gcirc_local = TgLginv(h);
		
	% Local drag, based on unit longitudinal drag, lateral according to the ratio, no local
	% torsional drag, multiplied by drag coefficient and local scaled
	% differential length
	gcirc_local_to_F_local = ...
        [-1     0       0;
        0   -drag_ratio 0;
        0       0       0]*c*lambda;
	
    % Transfer force to midpoint-tangent frame by transpose of the
    % adjoint-inverse action
    F_local_to_F_midpoint = transpose(Adjinv(h));
	
	% shape component of pfaffian
    % (map from shape velocities to midpoint-tangent forces)
	omega2 = F_local_to_F_midpoint ...
        * gcirc_local_to_F_local ...
        * gdot_to_gcirc_local ...
        * J; % J is rdot to gdot mapping
	
	% system body velocity component of pfaffian
	omega1 = F_local_to_F_midpoint ...
		* gcirc_local_to_F_local ...
        * gdot_to_gcirc_local ...
        * TeRg(h); % Mapping from gdot of system to gdot of outboard point
	
    % Combine contributions to Pfaffian from shape and position motion
	dOmega = [omega1 omega2];
    
    % Turn Pfaffian into column vector for ODE45
	dOmega = dOmega(:);

end