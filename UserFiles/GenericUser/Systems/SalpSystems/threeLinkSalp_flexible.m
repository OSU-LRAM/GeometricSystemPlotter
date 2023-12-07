function salp = threeLinkSalp_flexible(inputmode)

    %%%%%%
    %% Body geometry
    swimmerName = 'Three-Link Salp Continuous Curvature';

    switch inputmode
        case 'getName'
            salp = swimmerName;
        case 'getSwimmerSystem'

            s.name = swimmerName;
            %Chain of straight links
            s.geometry.type = 'general curvature';
            s.geometry.length = 1;
            s.geometry.baseframe = 'com-mean';
            s.geometry.function = @curv_PiecewiseCC_3Sections;
            s.geometry.nShapes = 3;
        
            %Set configuration space
            s.conf_space = LieGroups.SE2;
        
            %Set initial shape of salp for simulation
            s.geometry.initialShape = zeros(1,s.geometry.nShapes);
        
            %Set initial pose of salp for simulation
            s.geometry.initialPose = zeros(1,3);
        
            %% Thruster geometry
        
            %Specify whether arclength refers to portion of individual link or
            %portion of whole salp
            s.geometry.thrusters.arcParametrization = 'wholeChain';
            %Jet location represented as arc length on whole salp
            %-.5 is tip of tail, 0.5 is tip of head
            s.geometry.thrusters.thrusterArcLengths = [1/6:1/3:5/6]-.5; 
        
            thrusterMag = pi/6;
            thrusterSign = 2*mod(1:numel(s.geometry.thrusters.thrusterArcLengths),2)-1;
            s.geometry.thrusters.angles = thrusterMag*thrusterSign;
        
            s.geometry.thrusters.nThrusters = numel(s.geometry.thrusters.thrusterArcLengths);
        
            %% System physics
        
        
            %%%
            %%%%%%
            % Define system 
            s.physics.systemType = 'drag';
            s.physics.drag_ratio = 3;
            s.physics.drag_coefficient = 1;
            %s.physics.linear_drag_coefficient = 1;
            
            %Define system elasticity
            s.physics.jointStiffness = .005;
            %Top bit is elasticity of swimmer body wrt world frame: usually zero
            %Bottom bit is a matrix of spring connections.  For springs that do not
            %connect different shape modes, this is a diagonal matrix of spring
            %stiffnesses
            s.physics.K = [zeros(3,s.geometry.nShapes);diag(s.physics.jointStiffness*ones(1,s.geometry.nShapes))];
        
            %Functional distribution matrix
        
            s.ForceConnection = @(Shape, Force) salpForceToVelocity( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        Shape,...                                % Salp Shape
                        Force);                                  % Thrust Forces
        
            %% System visualization
        
            %shape space tic locations
            s.visualization.shape_tics = [-1 0 1];
        
            %%%% Video setup
            s.visualization.FrameRate = 30;
            % Buffer around salp centroid to plot as fraction of body length
            s.visualization.VidHalfWidth = 2/3;
        
            %%%% Grid settings
            s.visualization.showGrid = true;
            s.visualization.gridSpacing = 0.5;
        
        
            %%%% Animation colors
            s.visualization.colors.grid_color = [.8,.8,.8];  % Grey
            s.visualization.colors.link_color = [0.0745,0.6235,1.0000];  % Light Blue
            s.visualization.colors.thruster_color = [0.4941,0.1843,0.5569];  % Purple
            s.visualization.colors.in_jet_color = [1.0000,1.0000,0.0667];  % Yellow
            s.visualization.colors.out_jet_color = [1.0000,0.4118,0.1608];  % Orange
        
            %%%%
            %Save the system properties
            salp = s;
    end
end

