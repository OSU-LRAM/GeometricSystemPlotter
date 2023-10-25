function output = sysf_constantcurv3_lowRe(input_mode,pathnames)
% System file for a low Reynolds

    % Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Viscous swimmer: Constantcurve3'; % Display name

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
            s.geometry.function = {@(s)constant_curvature_3_1(s);@(s)constant_curvature_3_2(s);@(s)constant_curvature_3_3(s)};

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
            s.A = @(a1,a2,a3) LowRE_local_connection(s.geometry,s.physics,[a1;a2;a3]);
            
            % Locomotion model derived from viscous drag forces reacting to
            % local velocities of elements on the body
            s.metric = @(a1,a2,a3) LowRE_dissipation_metric(s.geometry,s.physics,[a1;a2;a3]);
            %%%


            %%%
            % Processing details

            %Range over which to evaluate connection
            s.grid_range = [-1,1,-1,1,-1,1]*5;

            %densities for various operations
            den = 6;
            s.density.vector = [den den den]; %density to display vector field
            s.density.scalar = [den den den]; %density to display scalar functions
            s.density.eval = [den den den];   %density for function evaluations
            s.density.metric_eval = [den den den];
            s.density.finite_element=den;

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