function [A, M_r] = tapeMeasure_connectionAndMetric(geometry,physics,shapeValues)

    %Get core values that define the shape of the fins
    C = geometry.tapeLength;
    alpha = shapeValues(1);
    L = shapeValues(2);
    
    %Get lengths of 2 fin components on either side of buckle point
    l1 = (C^2-L^2)/(2*C-2*L*cos(alpha));
    l2 = C - l1;
    
    %[X,Y] position of midpoint of links 1 and 2
    P1_top = [l1/2*cos(alpha),l1/2*sin(alpha)];
    P2_top = [l1/2*cos(alpha)+L/2,l1/2*sin(alpha)];
    
    %Symmetry gives coordinates of bottom midpoints
    P1_bottom = P1_top;
    P1_bottom(2) = -P1_top(2);
    P2_bottom = P2_top;
    P2_bottom(2) = -P2_top(2);
    
    %Take derivatives of x coordinate of link 1 wrt shape variables
    dP1tx_dalpha = -sin(alpha)*(C^2-L^2)/(4*(C-L*cos(alpha))) - L*cos(alpha)*sin(alpha)*(C^2-L^2)/(2*C-2*L*cos(alpha))^2;
    dP1tx_dL = cos(alpha)^2*(C^2-L^2)/(2*C-2*L*cos(alpha))^2 - L*cos(alpha)/(2*(C-L*cos(alpha)));
    
    %Derivative of y coordinate of link 1 wrt shape variables
    dP1ty_dalpha = cos(alpha)*(C^2-L^2)/(4*(C-L*cos(alpha))) - L*sin(alpha)^2*(C^2-L^2)/(2*C-2*L*cos(alpha))^2;
    dP1ty_dL = cos(alpha)*sin(alpha)*(C^2-L^2)/(2*C-2*L*cos(alpha))^2 - L*sin(alpha)/(2*(C-L*cos(alpha)));
    
    %Get beta (the angle between link2 and the controlled link) and
    %numerically find its derivative wrt shape variables
    beta = getBeta(alpha,C,L);
    step = 1e-5;
    dBeta_dalpha = (getBeta(alpha+step,C,L)-getBeta(alpha-step,C,L))/(2*step);
    dBeta_dL = (getBeta(alpha,C,L+step)-getBeta(alpha,C,L-step))/(2*step);
    
    %Start building jacobians and link mass matrices:
    
    %First is a central mass representing bouyancy platform and motors
    J0 = [eye(3),zeros(3,2)];
    m0 = diag([physics.headMass,physics.headMass,physics.headRotationalInertia]);
    
    %Second are the jacobians of the top and bottom first link, which are
    %symmetric about the controlled link
    J1t_body = [1 0 -P1_top(2) dP1tx_dalpha dP1tx_dL;...
        0 1 P1_top(1) dP1ty_dalpha dP1ty_dL;...
        0 0 1 1 0];
    J1b_body = [1 0 -P1_bottom(2) dP1tx_dalpha dP1tx_dL;...
        0 1 P1_bottom(1) -dP1ty_dalpha -dP1ty_dL;...
        0 0 1 -1 0];
    %X,Y mass depends on link length
    m1 = diag([physics.massRate,physics.massRate + physics.addedMassRate,0])*l1;
    %Rotational mass is link length to the third and fourth power for
    %inertia and added mass inertia, respectively
    m1(3,3) = physics.addedRotationalInertiaRate*(l1/2)^4 + physics.massRotationalInertiaRate*(l1/2)^3;
    
    %Finally we find the jacobians of the top and bottom second link, which
    %are symmetric about the controlled link
    J2t_body = [1 0 -P2_top(2) dP1tx_dalpha dP1tx_dL+1/2;...
        0 1 P2_top(1) dP1ty_dalpha dP1ty_dL;...
        0 0 1 -dBeta_dalpha -dBeta_dL];
    J2b_body = [1 0 -P2_bottom(2) dP1tx_dalpha dP1tx_dL+1/2;...
        0 1 P2_bottom(1) -dP1ty_dalpha -dP1ty_dL;...
        0 0 1 dBeta_dalpha dBeta_dL];
    %Find second link masses similar to the way we calculated first link
    %masses above
    m2 = diag([physics.massRate,physics.massRate + physics.addedMassRate,0])*l2;
    m2(3,3) = physics.addedRotationalInertiaRate*(l2/2)^4 + physics.massRotationalInertiaRate*(l2/2)^3;

    %These intermediate jacobians translate [x,y,theta] body velocities to
    %[x,y,theta] link velocities, which is the coordinate we need for our
    %geometry
    bodyRotation_1top = [cos(alpha),sin(alpha),0; -sin(alpha) cos(alpha) 0; 0 0 1];
    bodyRotation_1bottom = bodyRotation_1top';
    bodyRotation_2top = [cos(beta) -sin(beta) 0; sin(beta) cos(beta) 0; 0 0 1];
    bodyRotation_2bottom = bodyRotation_2top';
    
    %Convert each jacobian from body frame to link frame
    J1t = bodyRotation_1top*J1t_body;
    J1b = bodyRotation_1bottom*J1b_body;
    J2t = bodyRotation_2top*J2t_body;
    J2b = bodyRotation_2bottom*J2b_body;
    
    %Pull back each individual mass through the jacobians and sum to find
    %the total mass matrix
    M_full = J0'*m0*J0 + J1t'*m1*J1t + J1b'*m1*J1b + J2t'*m2*J2t + J2b'*m2*J2b;
    
    % Pfaffian is first three rows of M_full
    omega = M_full(1:3,:);
    % Build the local connection
    A = omega(:,1:3)\omega(:,4:end);
    A = real(A);
    % And fold down the mass matrix into the shape space
    M_r = [-A.' eye(size(A,2))] * M_full * [-A; eye(size(A,2))];
    M_r = real(M_r);
end

function beta = getBeta(alpha,C,L)

    l1 = (C^2-L^2)/(2*C-2*L*cos(alpha));
    l2 = C - l1;
    
    beta = acos((l2^2 + L^2 - l1^2)/(2*L*l2));

end