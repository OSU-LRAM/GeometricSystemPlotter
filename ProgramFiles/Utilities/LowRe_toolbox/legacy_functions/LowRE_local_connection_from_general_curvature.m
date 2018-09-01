function [A, h, J,Omega] = LowRE_local_connection_from_general_curvature(curvdef,cparams,L,c,drag_ratio)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
% curvdef: function mapping shape variables and position along body to curvature of body 
% cparams: value of shape variables
% L: total length of swimmer
% c: drag per unit length
% drag_ratio: ratio of lateral to longitudinal drag

% This function is now a wrapper on the general LowRE_local_connection
% function, maintained for legacy purposes.
geometry.type = 'general curvature';
geometry.function = curvdef;
geometry.length = L;

physics.drag_coefficient = c;
physics.drag_ratio = drag_ratio;


[A, h, J,Omega] = LowRE_local_connection(geometry,physics,cparams);

%warning('The function LowRE_local_connection_from_general_curvature for generating the local connection has been deprecated. For better maintainability and compatability with future sysplotter features, please update your sysf_ file to include a geometry structure, and call the generic LowRE_local_connection function with this structure.')


% %%%%%%%
% % Code that used to be run is below
% 
% % Specified integration limits
% int_limit = L*[-0.5 0.5];
% 
% % First, get the backbone locus and Jacobian functions
% [h, J] = backbone_from_general_curvature(curvdef,cparams,L);
% 
% % Now integrate to get the pfaffian
% Omega_sol = ode45( @(s,Omega) LowRE_Pfaffian_infinitesimal(s,h(s),J(s),c,drag_ratio),int_limit,zeros(3,3+length(cparams)));
% %Omega_sol = ode_multistart(@ode45,@(s,Omega) LowRE_Pfaffian_infinitesimal(s,h(s),J(s),c,drag_ratio),int_limit,int_limit(1),zeros(3,3+length(cparams)));
% 
% Omega = reshape(deval(Omega_sol,int_limit(end)),3,[]);
% %Omega = reshape(Omega_sol(int_limit(end)),3,[]);
% 
% A = Omega(:,1:3)\Omega(:,4:end);

end
