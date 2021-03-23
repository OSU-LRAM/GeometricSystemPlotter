function M_alpha = mass_matrix(geometry,physics,jointangles,s)
% Calculates the mass matrix necessary for optimizing the gaits of
% inertia-dominated swimmers.
% 
% Inputs:
%   geometry: Structure containing information about geometry of chain. See
%       Inertial_connection_discrete for information about required fields.
%
%   physics: Structure containing information about the system's physics.
%       See Inertial_connection_discrete for information about required 
%       fields.
%
%   jointangles: Array containing values of shape positions at which the
%       pulled-back partial mass matrix should be evaluated
%
% Outputs for system with k shape variables:
%   M_alpha: Pulled-back mass matrix, in terms of shape variables only. Of
%   size (k by k).

% Obtain the general mass matrix and local connection
[A, h, J, J_full, omega, M_full, local_inertias] = Inertial_connection_discrete(geometry,physics,jointangles);
% Pull back the mass matrix to be in terms of the shape variables only
M_alpha = [-A' eye(size(A,2))]*M_full*[-A; eye(size(A,2))];
end