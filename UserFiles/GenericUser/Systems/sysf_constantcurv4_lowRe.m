function output = sysf_constantcurv4_lowRe(input_mode,pathnames)
% System file for a low Reynolds

    % Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Viscous swimmer: Constantcurve4'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/ContinuousBackbone/',...
                'Physics/LowReynoldsRFT/'});

		case 'initialize'

            %%%%%%%%%
            % System definition

            %%%
            % Define the geometry of the system

            % this system's shape is defined by a set of basis functions
            % whose weighted sum is curvature as a function of distance
            % along the backbone
            s.geometry.type = 'curvature basis';                

            % The specific basis functions are those for a
            % piecewise-constant curvature system
            s.geometry.function = {@(s)constant_curvature_4_1(s);@(s)constant_curvature_4_2(s);@(s)constant_curvature_4_3(s);@(s)constant_curvature_4_4(s)};

            % Total length of the swimmer, in real units
            s.geometry.length = 1;
            
            % base the system off of its center frame
            s.geometry.baseframe = 'center';

            %%%
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1 -0.5 0 0.5 1]*5;
            
            %%%

            %%%
            % Define the physics of the system

            % This system is treated as moving in a viscous fluid. Key elements we
            % need for this model are a drag coefficient and a drag ratio
            s.physics.drag_coefficient = 1;                   % multiplier from longitudinal velocity to drag force. Changing this scales the dissipation matrix
            s.physics.drag_ratio = 2;                         % ratio of lateral:longitudinal drag

            % Locomotion model derived from viscous drag forces reacting to
            % local velocities of elements on the body
            s.A = @(a1,a2,a3,a4) LowRE_local_connection(s.geometry,s.physics,[a1;a2;a3;a4]);
            
            % Locomotion model derived from viscous drag forces reacting to
            % local velocities of elements on the body
            s.metric = @(a1,a2,a3,a4) LowRE_dissipation_metric(s.geometry,s.physics,[a1;a2;a3;a4]);
            %%%


            %%%
            % Processing details

            %Range over which to evaluate connection
            s.grid_range = [-1,1,-1,1,-1,1,-1,1]*6;

            %densities for various operations
            s.density.vector = [13 13 13 13]; %density to display vector field
            s.density.scalar = [13 13 13 13]; %density to display scalar functions
            s.density.eval = [13 13 13 13];   %density for function evaluations
            s.density.metric_eval = [1 1 1 1]*13;
            s.density.finite_element =7;
            %shape space tic locations
            s.tic_locs.x = [-1 0 1]*6;
            s.tic_locs.y = [-1 0 1]*6;
            %%%%
            
            % Set system type variable for gait optimization
            s.system_type = 'drag';
    
    
			%Save the system properties
			output = s;


	end
end