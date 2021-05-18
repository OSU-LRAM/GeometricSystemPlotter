function [A, h, J, J_full, omega] = LowRE_connection_discrete(geometry,physics,jointangles)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields required for this function are:
%
%       linklengths: A vector of link lengths defining a kinematic chain
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
%           middle link, or at the joint between the two middle links and at
%           the mean of their orientation. Options for this field are
%           specified in the documentation of N_link_chain.
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%          is not provided or is entered as an empty matrix, then the links
%          will not be scaled.
%
%   physics: Structure containing information about the system's physics.
%       Fields are:
% 
%       drag_ratio: Ratio of lateral to longitudinal drag
%
%       drag_coefficient: Drag per unit length for longitudinal direction
%
%   jointangles: Angles of the joints between the links, defining the
%       chain's current shape.
%
% Outputs:
%
%   A: The "local connection" matrix, or "Locomotion Jacobian". This matrix
%       maps joint angular velocities ("shape velocities") to body
%       velocities of the system's base frame as 
%
%            g_b = -A * alphadot
%
%       with g_b the body velocity and alphadot the joint angular velocity.
%
%       (Note the negative sign in this equation; sysplotter code uses the
%       classical formulation of the geometric mechanics equations, which
%       includes this negative sign.
%
%   [h,J,J_full]: location and Jacobians for the links on the chain. These
%       are passed through from N_link_chain; see the documentation on that
%       function for more details
%
%   omega: The "Pfaffian constraint matrix" from which the local connection
%       A is calculated. This matrix is a linear map from system body and
%       shape velocities to net external forces acting on the base frame,
%       which must be zero for all achievable motions of the system


    %%%%
    % First, get the positions of the links in the chain and their
    % Jacobians with respect to the system parameters
    if isfield(geometry,'subchains')
        [h, J, J_full] = branched_chain(geometry,jointangles);
    else
        [h, J, J_full] = N_link_chain(geometry,jointangles);
    end

    %%%%%%%%
    % We are modeling low Reynolds number physics as being resistive
    % viscous force, with a quasistatic equilibrium condition that inertial
    % forces are negligable, so that the viscous drag on the whole system
    % is zero at any given time. (The drag is anisotropic, so the system
    % can actually be moving even though there is no net force acting on
    % it).
    %
    % Because viscous drag is linear, we can build a linear map from the
    % body velocity of a given link to the force and moment acting on the
    % link.
    % 
    % Further, because the velocity kinematics from the system's body and
    % shape velocities to the body velocity of each link are linearly
    % encoded by the Jacobians of the links, we can use these Jacobians to
    % pull the drag operators on the links back onto the space of system
    % body and shape velocities, giving us a linear map from the system's
    % body and shape velocity to the net forces acting on the system
    %    
    % To calculate this linear mapping, we first pre-allocate storage for
    % the metric contributions. Each contribution is 3 x m (rows for
    % x,y,theta motion, and one column per joint), and there is one
    % contribution per link. This structure is of the same dimensions as
    % J_full, so we use it as a template.
    link_force_maps = J_full;
    
    % Now iterate over each link, calculating the map from system body and
    % shape velocities to forces acting on the body
    for idx = 1:numel(link_force_maps)
        
        link_force_maps{idx}= LowRE_body_drag_link(h.pos(idx,:),...            % Position of this link relative to the base frame
                                                    J_full{idx},...             % Jacobian from body velocity of base link and shape velocity to body velocity of this link
                                                    h.lengths(idx),...          % Length of this link
                                                    physics.drag_ratio,...      % Ratio of lateral to longitudinal drag
                                                    physics.drag_coefficient);  % Bulk drag coefficient
  
    end

    % Sum the force-maps for each link to find the total map from system
    % body and shape velocities to force actign on the body
    omega = sum(cat(3,link_force_maps{:}),3);

    %%%%%%%%
    % Because of the quasistatic equilibrium condtion, if there are no
    % outside forces acting on the swimmer, the viscous forces must be zero
    % during all of its motions, and thus allowable motions for the system
    % are in the nullspace of the omega matrix (i.e. omega acts as a
    % "Pfaffian constraint"),
    %
    % [0] = omega * [g_b; alphadot]
    %
    % For our locomotion studies, we will be specifying shape velocities and
    % mapping them to the body velocities that satisfy the constraint
    % equation. To construct this map, we split omega into the blocks that
    % act on g_b and alphadot,
    %
    % [0] = [omega_g omega_alpha] * [g_b; alphadot],
    %
    % pull the body velocity terms to the left,
    %
    % -(omega_g * g_b) = omega_alpha * alphadot,
    %
    % and then invert the omega_g term,
    %
    % g_b = - ( inv(omega_g) * omega_alpha ) * alphadot.
    %
    % Grouping the two omega terms then gives  the local connection as
    % 
    % g_b = - A * alphadot
    %
    % (Where the notational choice to leave the negative sign out of the
    % grouping is an historical artifact from the history of geometric
    % mechanics)
    
    % For this calculation, omega_g is the first three columns of omega,
    % and omega_alpha is the remaining columns
    A = omega(:,1:3)\omega(:,4:end);


end


function omega = LowRE_body_drag_link(h,J_full,L,drag_ratio,c)
% Calculate the matrix that maps from system body and shape velocities to
% forces acting on the base frame of the system

		
	%%%%%%%
    % Forces acting on the system are applied as forces acting on the links.
    %
    % Body forces on a link are mapped to body forces on the base frame by
    % the transpose (or "dual") of the Adjoint-inverse action of the
    % position of the link relative to the base frame
    
    % F_local is body force acting on the link
    % F_midpoint is body force acting on the baseframe
    F_local_to_F_baseframe = transpose(Adjinv(h));
    
    
    %%%%%%%%
    % In a viscous medium, the forces acting on a link can be represented
    % as a matrix that multiplies the body velocity of the link to produce
    % body forces acting on the link.
    %
    % Under a simple resistive-force model, this drag matrix is diagonal,
    % with different coefficients for longitudinal, lateral, and rotational
    % velocity.
    %
    % For our model:
    %
    %   Longitudinal drag is proportional to longitudinal velocity and link
    %       length
    %
    %   Lateral drag is proportional to lateral velocity and link length,
    %       with a different drag coefficient, specified by its ratio to
    %       the longitudinal drag.
    %
    %   Rotational drag is proportional to rotational velocity, with a
    %       coefficient based on integrating lateral drag along the length
    %       of the link as it rotates: 
    %
    %           Moment = thetadot * int_{-L/2}^{L/2} (drag_ratio * s^2) ds
    %                  = drag_ratio/12 * L^3
    %
    %   c is the absolute drag coefficient for the body in the fluid
    
    % Drag matrix with terms as describe above.
    % gcirc_local is the body velocity of the link
    % F_local is the body force acting on the link
    gcirc_local_to_F_local = ...
        [-L      0               0;
        0    -drag_ratio*L       0;
        0        0           -drag_ratio/12*L^3]*c;
    
    %%%%%%%%%%
    % Mapping from system body and shape velocity to body velocity of the
    % link is given by J_full
    %
    % Linear map from body and shape velocity to force on body frame is
    % therefore
    % product of the two matrices calculated above and J_full
    omega = F_local_to_F_baseframe * gcirc_local_to_F_local * J_full;

end