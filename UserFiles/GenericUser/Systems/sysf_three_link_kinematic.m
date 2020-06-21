function output = sysf_three_link_kinematic(input_mode,pathnames)


	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end


	
	switch input_mode

		case 'name'

			output = 'Three-link Kinematic Snake'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/Nonholonomic/'});

		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = [1 1 1];
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            s.geometry.constraint_list = struct('link_number',{1;2;3},...               % Links on which the wheels are placed
                                     'constraint_direction',{[0 1 0];[0 1 0];[0 1 0]}); % Orientations of the wheels
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1  0  1];
            
            % Tell sysplotter to draw this system with wheels on the links
            s.visual.drawing_function = @wheeled_chain;

            
            %%%%%

			%Functional representation of local connection, split into a
			%numerator and denominator because the connection has
			%singularities

			s.A_num = @(a1,a2) Nonholonomic_connection_discrete(s.geometry,[a1,a2],'num');
			s.A_den = @(a1,a2) Nonholonomic_connection_discrete(s.geometry,[a1,a2],'den');


			%%%
			%Processing details

			%Mark that system has a singularity that needs to be accounted for
			s.singularity = 1;

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*pi/2;

			%densities for various operations
			s.density.vector = [10 10]; %density to display vector field
			s.density.scalar = [30 30]; %density to display scalar functions
			s.density.eval = [30 30];   %density for function evaluations
            s.density.finite_element=30;

			%%%
			%Display parameters

			%shape space tic locations
			s.tic_locs.x = [-1 0 1];
			s.tic_locs.y = [-1 0 1];

            % Set system type variable for gait optimization
            s.system_type = 'inertia';

			%%%%
			%Save the system properties
			output = s;



	end

end