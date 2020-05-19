addpath('Physics/Inertial/Discrete_links');
addpath('Physics/Inertial');
addpath('Utilities/Kinematics');
syms a1 a2;

link_shape_params = struct('aspect_ratio',0.125);
geometry.link_shape_parameters = repmat({link_shape_params},3,1);
link_shape = struct('a',0,'b',0);
geometry.link_shape = repmat({link_shape},3,1);
geometry.type = 'n-link chain';
geometry.linklengths = [1,1,1];
geometry.function = '';
physics = struct('drag_ratio','','drag_coefficient','','fluid_density',1);
shapeparams = [a1,a2];
[A, h, J,J_full,Omega] = Inertial_local_connection(geometry,physics,shapeparams);

% Get the partial derivative of the Jacobian
dJdq = mobile_jacobian_derivative(J_full);

% Now iterate over each link, calculating the map from system body and
% shape velocities to forces acting on the body
for idx = 1:numel(J_full)

    [link_inertias{idx},local_inertias{idx}] = Inertia_link(h.pos(idx,:),...            % Position of this link relative to the base frame
                                                J_full{idx},...             % Jacobian from body velocity of base link and shape velocity to body velocity of this link
                                                h.lengths(idx),...          % Length of this link
                                                geometry.link_shape{idx},...         % Shape type of this link
                                                geometry.link_shape_parameters{idx},...  % Shape parametes for this link
                                                physics.fluid_density);      % Fluid density relative to link

end
    
% Get the partial derivative of the Mass matrix
dMdq = partial_mass_matrix(J_full,dJdq,local_inertias,'mobile');