function output = sysf_two_mode_serpenoid(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'HighRe Serpenoid, Two-Mode'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/Inertial/'});
            
		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = ones(1,20);
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            s.geometry.link_shape = {};
            st = struct('aspect_ratio',0.1);
            for i = 1:20
                s.geometry.link_shape{i} = 'ellipse';
                s.geometry.link_shape_parameters{i} = st;
            end
            modeSetup = [];
            for i = 1:19
                modeSetup(i,:) = [cos(i/20*2*pi),sin(i/20*2*pi)];
            end
            s.geometry.modes = modeSetup;
            
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1  0  1]*.5;
            
            %%%
            %%%%%%
            % Define system physics
            s.physics.fluid_density = 1;
         
            %Functional Local connection and dissipation metric


            s.A = @(alpha1,alpha2) Inertial_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2]);                        % Joint angles

            s.metric = @(alpha1,alpha2) Inertial_energy_metric(s.geometry,s.physics,[alpha1,alpha2]);
            %s.metric = @(alpha1,alpha2) eye(2);
                    
			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*.5;

			%densities for various operations
			s.density.vector = [21 21]; %density to display vector field
			s.density.scalar = [51 51]; %density to display scalar functions
			s.density.eval = [31 31];   %density for function evaluations
            s.density.metric_eval = [11 11]; %density for metric evaluation
            s.density.mass_eval = [31 31]; % density for mass matrix evaluation
            s.density.coriolis_eval = [31 31];
            s.density.finite_element=31;


			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*.5;
			s.tic_locs.y = [-1 0 1]*.5;


            % System type flag
            s.system_type = 'inertia';
			%%%%
			%Save the system properties
			output = s;


	end

end

