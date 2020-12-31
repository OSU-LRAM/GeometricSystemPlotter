function output = sysf_four_link_floating_ellipse(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Inertial floating: 4-link ellipses'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/Inertial/'});
            
		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = [1 1 1 1];
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            s.geometry.link_shape = {'ellipse','ellipse','ellipse','ellipse'};
                st = struct('aspect_ratio',0.01);
            s.geometry.link_shape_parameters = {st,st,st,st};
            
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1  0  1]*1;
            
            %%%
            %%%%%%
            % Define system physics
            s.physics.fluid_density = 0;
           
 
            %Functional Local connection and dissipation metric

            s.A = @(alpha1,alpha2,alpha3) Inertial_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2,alpha3]);                        % Joint angles
            
            s.metric = @(alpha1,alpha2,alpha3) Inertial_energy_metric(s.geometry,s.physics,[alpha1,alpha2,alpha3]);
            %s.metric = @(alpha1,alpha2,alpha3) eye(3);
             
            % TODO: These should probably be calculated as part of a larger
            % wrapping function that's meant to return M and C matrices for
            % a set of points
%             s.dJdq = @(alpha1,alpha2) mobile_jacobian_derivative(s.J_full);
%             s.dMdq = @(alpha1,alpha2) partial_mass_matrix(s.J,s.dJdq,local_inertias,'mobile');
%             s.M_alpha = @(alpha1,alpha2) mass_matrix(s.geometry,s.physics,[alpha1,alpha2]);
%             s.dM_alphadalpha = @(alpha1,alpha2,A_eval,A_grid) shape_partial_mass(s.geometry,s.physics,[alpha1,alpha2],A_eval,A_grid);
%  

			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1,-1,1]*2*pi/3;

			%densities for various operations
            sampleDensity = 12;
            sd = sampleDensity;
			s.density.vector = [sd sd sd]; %density to display vector field
			s.density.scalar = [sd sd sd]; %density to display scalar functions
			s.density.eval = [sd sd sd];   %density for function evaluations
            s.density.metric_eval = [sd sd sd]; %density for metric evaluation
            s.density.mass_eval = [sd sd sd]; % density for mass matrix evaluation
            s.density.coriolis_eval = [sd sd sd];
            s.density.metric_display = [sd sd sd];
            s.density.finite_element=sd;


			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*1;
			s.tic_locs.y = [-1 0 1]*1;

            % System type flag
            s.system_type = 'inertia';
			%%%%
			%Save the system properties
			output = s;


	end

end

