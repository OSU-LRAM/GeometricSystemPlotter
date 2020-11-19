function output = sysf_three_link_lowRe(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Viscous Swimmer: 3-link'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/LowReynoldsRFT/'});
            
		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = [1 1 1];
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            s.geometry.link_shape = {'ellipse','ellipse','ellipse'};
                st = struct('aspect_ratio',0.1);
            s.geometry.link_shape_parameters = {st,st,st};
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1  0  1];

            %%%
            %%%%%%
            % Define system physics
            s.physics.drag_ratio = 2;
            s.physics.drag_coefficient = 1;
            s.physics.fluid_density = 1;
           
            %Functional Local connection and dissipation metric

            s.A = @(alpha1,alpha2) LowRE_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2]);                        % Joint angles
            
            s.metric = @(alpha1,alpha2) LowRE_dissipation_metric(...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2]);                        % Joint angles
     
			%%%
            
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*2.5;

			%densities for various operations
			s.density.vector = [21 21 ]; %density to display vector field
			s.density.scalar = [51 51 ]; %density to display scalar functions
			s.density.eval = [31 31 ];   %density for function evaluations
            s.density.metric_eval = [11 11]; %density for metric evaluation
            s.density.finite_element=31;
            s.density.coriolis_eval = [31 31];
            s.density.mass_eval = [31 31]; % density for mass matrix evaluation

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*1;
			s.tic_locs.y = [-1 0 1]*1;

            % Set system type variable for gait optimization
            s.system_type = 'drag';
            
			%%%%
			%Save the system properties
			output = s;


	end

end

