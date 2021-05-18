function output = sysf_three_mode_CC(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'HighRe CC, 3 Mode'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/Inertial/'});
            
		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = ones(1,19);
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            s.geometry.link_shape = {};
            st = struct('aspect_ratio',0.1);
            for i = 1:19
                s.geometry.link_shape{i} = 'ellipse';
                s.geometry.link_shape_parameters{i} = st;
            end
            modeSetup = [];
            for i = 1:18
                if i <= 6
                    modeSetup(i,:) = [1/40,0,0];
                elseif i <= 12
                    modeSetup(i,:) = [0,1/40,0];
                else
                    modeSetup(i,:) = [0,0,1/40];
            end
            s.geometry.modes = modeSetup;
            
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1  0  1]*15;
            
            %%%
            %%%%%%
            % Define system physics
            s.physics.fluid_density = 1;
            %s.physics.interaction = 'off';
 
            %Functional Local connection and dissipation metric


            s.A = @(alpha1,alpha2,alpha3) Inertial_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2,alpha3]);                        % Joint angles

            s.metric = @(alpha1,alpha2,alpha3) Inertial_energy_metric(s.geometry,s.physics,[alpha1,alpha2,alpha3]);
            %s.metric = @(alpha1,alpha2,alpha3) eye(3);
             
			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1,-1,1]*20;
            
            sampleDensity = 12;
            sd = sampleDensity;
			%densities for various operations
			s.density.vector = [sd sd sd]; %density to display vector field
			s.density.scalar = [sd sd sd]; %density to display scalar functions
			s.density.eval = [sd sd sd];   %density for function evaluations
            s.density.metric_eval = [sd sd sd]; %density for metric evaluation
            s.density.mass_eval = [sd sd sd]; % density for mass matrix evaluation
            s.density.coriolis_eval = [sd sd sd];
            s.density.finite_element=sd;


			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*15;
			s.tic_locs.y = [-1 0 1]*15;


            % System type flag
            s.system_type = 'inertia';
			%%%%
			%Save the system properties
			output = s;


	end

end

