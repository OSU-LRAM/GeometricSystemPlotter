%Converts forces applied at thrust locations on a salp at a given shape to 
%resultant body and shape velocities for a drag dominated swimmer
function [g_circ,shapeVel] = salpForceToVelocity_DragDominated(geometry,physics,shape,forces)

    switch geometry.type
        
        case 'n-link chain'
            %Get jacobians for thruster locations and link centers
            chainOptions.useLinkCenters = false;
            chainOptions.getJacobianDerivative = false;
            [~,~,J_full,~,~,~,J_full_linkCenters] = N_link_chain(geometry,shape,chainOptions);
            %Get distribution matrix mapping link forces to body/shape forces
            B = getDistribution_discrete(J_full);
            %Get drag matrix from link centers
            D = getDragMatrix(J_full_linkCenters,geometry,physics);
            %Get stiffness matrix from system physics
            K = physics.K;

            %Solve force equation relating thrust and spring forces to drag
            %forces to get resulting body motion
            allVelocities = inv(D)*(B*forces(:) - K*shape(:));
            %Break these results into body velocities
            g_circ = allVelocities(1:3);
            %And shape velocities
            shapeVel = allVelocities(4:end);

        case 'general curvature'
            [~,~,D,J_full] = LowRE_connection_and_metric_continuous(geometry,physics,shape);

            J_full_cell = cell(1,geometry.thrusters.nThrusters);
            for i = 1:geometry.thrusters.nThrusters
                J_full_cell{i} = J_full(geometry.thrusters.thrusterArcLengths(i));
            end

            B = getDistribution_discrete(J_full_cell);
            K = physics.K;

            allVelocities = inv(D)*(B*forces(:) - K*shape(:));
            g_circ = allVelocities(1:3);
            shapeVel = allVelocities(4:end);
        otherwise
            error(['Whatever salp geometry you''re trying to use hasn''t been implemented yet or is named incorrectly.']);
    end

end
         

%Gets a salp drag matrix using the Jacobian for each link and definitions
%for salp geometry and salp physics
function D = getDragMatrix(J_lc,geometry,physics)

    %Drag coefficient
    cd = physics.drag_coefficient;
    %Link lengths, normalized to total salp length
    linklengths = geometry.linklengths/(sum(geometry.linklengths))*geometry.length;
    %Ratio of lateral drag to forward drag
    dragRatio = physics.drag_ratio;

    %Initialize drag matrix to zeros of correct dimension
    D = zeros(size(J_lc{1},2));

    %Sum drag contribution from each link
    for i = 1:numel(J_lc)
        J = J_lc{i};
        L = linklengths(i);

        drag_matrix =     [L       0                0;...
                    0    dragRatio*L       0;...
                    0        0           dragRatio/12*L^3]*cd;...

        D = D + J'*drag_matrix*J;
    end

end