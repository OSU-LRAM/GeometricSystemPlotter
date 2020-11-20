function output = sysf_three_link_HighRe(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = '5-Link Flapper with Passive Joints'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/Inertial/'});
            
		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = [1 1 1 1 1];
            s.geometry.activelinks = [1 1 1 1 1];
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            s.geometry.link_shape = {'ellipse','ellipse','ellipse','ellipse','ellipse'};
                st = struct('aspect_ratio',0.1);
            s.geometry.link_shape_parameters = {st,st,st,st,st};
            s.geometry.modes = [1,0;0,1;0,1;1,0];
            s.geometry.linkSeparation = .1;
            
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
            s.physics.interaction = 'off';
           
 
            %Functional Local connection and dissipation metric

            s.A = @(alpha1,alpha2) Inertial_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2]);                        % Joint angles
            
            s.metric = @(alpha1,alpha2)eye(2);%@(alpha1,alpha2) LowRE_dissipation_metric(...
%                         s.geometry,...                           % Geometry of body
%                         s.physics,...                            % Physics properties
%                         [alpha1,alpha2]);                        % Joint angles

            % TODO: These should probably be calculated as part of a larger
            % wrapping function that's meant to return M and C matrices for
            % a set of points
%             s.dJdq = @(alpha1,alpha2) mobile_jacobian_derivative(s.J_full);
%             s.dMdq = @(alpha1,alpha2) partial_mass_matrix(s.J,s.dJdq,local_inertias,'mobile');
                    
			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1.2,1.2,-2,2];

			%densities for various operations
			s.density.vector = [21 21]; %density to display vector field
			s.density.scalar = [51 51]; %density to display scalar functions
			s.density.eval = [31 31];   %density for function evaluations
            s.density.metric_eval = [11 11]; %density for metric evaluation
            s.density.mass_eval = [31 31]; % density for mass matrix evaluation
            s.density.coriolis_eval = [31 31];
            s.density.finite_element=31;


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

