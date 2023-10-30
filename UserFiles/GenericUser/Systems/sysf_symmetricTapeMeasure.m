function output = sysf_symmetricTapeMeasure(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Symmetric Tape Measure Swimmer: High-Re'; % Display name

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
            %.56 is real tape length in m
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
            
            s.physics.massRate = 0.001;
            s.physics.massRotationalInertiaRate = .01;
            
            s.physics.addedMassRate = 1;
            s.physics.addedRotationalInertiaRate = s.physics.addedMassRate/8;

            s.physics.headMass = .25;
            s.physics.headRotationalInertia = 1;
            
           
 
            %Functional Local connection and dissipation metric


            s.A = @(alpha1,alpha2) tapeMeasure_localConnection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2]);                        % Joint angles

            s.metric = @(alpha1,alpha2) tapeMeasure_metric(s.geometry,s.physics,[alpha1,alpha2]);
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
			s.grid_range = [0,pi/2,.24/.56,.45/.56];

			%densities for various operations
			s.density.vector = [21 21]; %density to display vector field
			s.density.scalar = [51 51]; %density to display scalar functions
			s.density.eval = [31 31];   %density for function evaluations

            %Bump up this value
            s.density.metric_eval = [31 31]; %density for metric evaluation
%            s.density.mass_eval = [31 31]; % density for mass matrix evaluation
            s.density.metric_display = [7 7];
%            s.density.coriolis_eval = [31 31];
            s.density.finite_element=31;


			%shape space tic locations
			s.tic_locs.x = [0,pi/4];
			s.tic_locs.y = [.24/.56,.45/.56];


            % System type flag
            s.system_type = 'inertia';
			%%%%
			%Save the system properties
			output = s;


	end

end

