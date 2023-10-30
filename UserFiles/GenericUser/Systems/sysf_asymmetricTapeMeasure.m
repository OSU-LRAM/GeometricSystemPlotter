function output = sysf_asymmetricTapeMeasure(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Asymmetric Tape Measure Swimmer: High-Re'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/Inertial/TapeMeasure'});
            
		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'looped chain';
            s.geometry.linklengths = [1 1 1 1];
            s.geometry.baseframe = 'center';
            s.geometry.tapeLength = 1;
            
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1  0  1];
            
            %%%
            %%%%%%
            % Define system physics
            s.physics.fluid_density = 1000;
            s.physics.addedMassRate = 0.2;
            s.physics.massRate = 0.001;
            s.physics.addedRotationalInertiaRate = 5.9;
            s.physics.massRotationalInertiaRate = 1/15;
            s.physics.headMass = .25;
            s.physics.headRotationalInertia = 1;
            
           
 
            %Functional Local connection and dissipation metric


            s.A = @(alpha1,alpha2,bodyLength) asymTapeMeasure_localConnection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2,bodyLength]);                        % Joint angles

            s.metric = @(alpha1,alpha2,bodyLength) asymTapeMeasure_metric(s.geometry,s.physics,[alpha1,alpha2,bodyLength]);
            %s.metric = @(alpha1,alpha2) eye(2);
             
%             % TODO: These should probably be calculated as part of a larger
%             % wrapping function that's meant to return M and C matrices for
%             % a set of points
% %             s.dJdq = @(alpha1,alpha2) mobile_jacobian_derivative(s.J_full);
% %             s.dMdq = @(alpha1,alpha2) partial_mass_matrix(s.J,s.dJdq,local_inertias,'mobile');
%             s.M_alpha = @(alpha1,alpha2) mass_matrix(s.geometry,s.physics,[alpha1,alpha2]);
%             s.dM_alphadalpha = @(alpha1,alpha2,A_eval,A_grid) shape_partial_mass(s.geometry,s.physics,[alpha1,alpha2],A_eval,A_grid);
                    
			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [0,pi/3,0,pi/3,24/56,45/56];

			%densities for various operations
            sample_density = 12;
            sd = sample_density;
			s.density.vector = [sd sd sd]; %density to display vector field
			s.density.scalar = [sd sd sd]; %density to display scalar functions
			s.density.eval = [sd sd sd];   %density for function evaluations
            %Bump up this value
            s.density.metric_eval = [sd sd sd]; %density for metric evaluation
%            s.density.mass_eval = [31 31]; % density for mass matrix evaluation
            s.density.metric_display = [7 7 7];
%            s.density.coriolis_eval = [31 31];
            s.density.finite_element=sd;


			%shape space tic locations
			s.tic_locs.x = [0,pi/3];
            s.tic_locs.y = [0,pi/3];
			s.tic_locs.z = [24/56,45/56];


            % System type flag
            s.system_type = 'inertia';
			%%%%
			%Save the system properties
			output = s;


	end

end

