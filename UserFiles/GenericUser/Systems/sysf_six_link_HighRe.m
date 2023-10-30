function output = sysf_six_link_HighRe(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'HighRe ideal swimmer: 6-link ellipses'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/Inertial/'});
            
		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = [1 1 1 1 1 1];
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            s.geometry.link_shape = {'ellipse','ellipse','ellipse','ellipse','ellipse','ellipse'};
                st = struct('aspect_ratio',0.1);
            s.geometry.link_shape_parameters = {st,st,st,st,st,st};
            %s.geometry.modes = [sqrt(2)/2,-sqrt(2)/2;sqrt(2)/2,sqrt(2)/2];
            
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1  0  1];
            
            %%%
            %%%%%%
            % Define system physics
            s.physics.fluid_density = 1;
           
 
            %Functional Local connection and dissipation metric


            s.A = @(alpha1,alpha2,alpha3,alpha4,alpha5) Inertial_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2,alpha3,alpha4,alpha5]);                        % Joint angles

            s.metric = @(alpha1,alpha2,alpha3,alpha4,alpha5) Inertial_energy_metric( ...
                        s.geometry,s.physics,[alpha1,alpha2,alpha3,alpha4,alpha5]);
            %s.metric = @(alpha1,alpha2,alpha3,alpha4) eye(4);
             
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
			s.grid_range = [-1,1,-1,1,-1,1,-1,1,-1,1]*1.5;

			%densities for various operations
            sample_density = 5;
            sd = sample_density;
			s.density.vector = [sd sd sd sd sd]; %density to display vector field
			s.density.scalar = [sd sd sd sd sd]; %density to display scalar functions
			s.density.eval = [sd sd sd sd sd];   %density for function evaluations
            %Bump up this value
            s.density.metric_eval = [sd sd sd sd sd]; %density for metric evaluation
%            s.density.mass_eval = [31 31]; % density for mass matrix evaluation
            s.density.metric_display = [sd sd sd sd sd];
%            s.density.coriolis_eval = [31 31];
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

