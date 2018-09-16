function [A, h, J, Omega] = LowRE_local_connection_from_curvature_bases(kappa_basis_input,r,L,c,drag_ratio)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
% kappa_basis_input: cell vector of handles to basis functions
% r: coefficients for basis functions
% L: total length of swimmer
% c: drag per unit length
% drag_ratio: ratio of lateral to longitudinal drag


% This function is now a wrapper on the general LowRE_local_connection
% function, maintained for legacy purposes.
geometry.type = 'curvature basis';
geometry.function = kappa_basis_input;
geometry.length = L;

physics.drag_coefficient = c;
physics.drag_ratio = drag_ratio;

[A, h, J,Omega] = LowRE_local_connection(geometry,physics,r);

%warning('The function LowRE_local_connection_from_curvature_bases for generating the local connection has been deprecated. For better maintainability and compatability with future sysplotter features, please update your sysf_ file to include a geometry structure, and call the generic LowRE_local_connection function with this structure.')

%%%%%%%
% Code that used to be run is below

% % Specified integration limits
% int_limit = L*[-0.5 0.5];
% 
% % First, get the backbone locus and Jacobian functions
% [h, J] = backbone_from_curvature_bases(kappa_basis_input,r,L);
% 
% % Now integrate to get the Pfaffian
% Omega_sol = ode45( @(s,Omega) LowRE_Pfaffian_infinitesimal(s,h(s),J(s),c,drag_ratio),int_limit,zeros(3,3+length(r)));
% 
% % Reshape the terms of the Pfaffian into a matrix of the correct dimension
% Omega = reshape(deval(Omega_sol,int_limit(end)),3,[]);
% 
% % Calculate the local connection by multiplying the inverse of the first
% % block in the Pfaffian by the second block
% A = Omega(:,1:3)\Omega(:,4:end);
% 
end


