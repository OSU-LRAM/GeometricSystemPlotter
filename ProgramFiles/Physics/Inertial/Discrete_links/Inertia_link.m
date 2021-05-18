function [Inertia_link_system,Inertia_link_local] = Inertia_link(h,J_full,L,link_shape_parameters,fluid_density)
% Calculate the matrix that maps from system body and shape velocities to
% forces acting on the base frame of the system

    % These are for an elliptical link. Can be generalized to other shapes

    % length and width of an elliptical link
    a = L/2;
    b = link_shape_parameters.aspect_ratio * a;
    
	mass_link = pi*a*b;
    
    rotational_inertia_link = mass_link*(a^2+b^2)/4;
    
    added_mass = fluid_density * diag([pi*b^2, pi*a^2, (a^2-b^2)^2]);
    
    Inertia_link_local = diag([mass_link,mass_link,rotational_inertia_link]) + added_mass;
    
    
    % Pullback of local inertia to body velocity/shape velocity coordinates
    
    Inertia_link_system = J_full' * Inertia_link_local * J_full;

end