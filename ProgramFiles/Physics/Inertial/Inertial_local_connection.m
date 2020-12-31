%function [A, M_a,J_full, local_inertias,M_full,h,J] = Inertial_local_connection(geometry,physics,shapeparams)
function [A, h, J,J_full,Omega,M_full] = Inertial_local_connection(geometry,physics,shapeparams)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
% geometry: structure defining system geometry
%      geometry.type: how the system geometry is defined 
%         (e.g., links, curvature basis, general curvature)
%      geometry.function: map from shape variables to local backbone
%         deformation (e.g., curvature or joint angles)
%      geometry.length: total length of swimmer
% physics: structure defining system physics
%      fluid_density: fluid density relative to body density
%      interactions: boolean as to whether the added mass interactions
%           between body panels should be considered
% cparams: value of shape variables




% Identify what kind of system is being calculated, and use this to specify how
% the local connection should be generated
switch geometry.type
    
    case {'curvature basis','curvature bases','general curvature'}
        physics_function = @Inertial_connection_continuous;
        
    case {'n-link chain','branched chain'}
        physics_function = @Inertial_connection_discrete;
        
end

% Call the physics function identified for the system
[A, h, J,Omega,J_full,M_full] = physics_function(geometry,physics,shapeparams);

end