function output = sysf_fishTail_HighRe_untapered(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

        
        
		case 'name'

			output = 'Flagella Tail Swimmer - Untapered'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/NLinkChain/',...
                'Physics/Inertial/'});
            
		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            nLinksTail = 20;
            nLinksHead = 10;
            nLinks = nLinksTail + nLinksHead;
            s.geometry.linklengths = [1/nLinks*ones(1,nLinks)];

            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            s.geometry.link_shape = {};
            st = struct('aspect_ratio',0.1);
            for i = 1:nLinks
                s.geometry.link_shape{i} = 'ellipse';
                s.geometry.link_shape_parameters{i} = st;
            end

            tailMode = [];
            nJointsTail = nLinksTail-1;
            C = (12*nJointsTail + 12)/(3*nJointsTail^2 + 2*nJointsTail);
            for i = 1:nJointsTail
                x = i/nLinksTail;
                tailMode(i,:) = [C*(x^2 - 1/3*x^3),0];
            end
            headMode = [0,1];
            rigidMode = zeros(nLinksHead-1,2);

            modeSetup = [tailMode;headMode;rigidMode];
            s.geometry.modes = modeSetup;
            
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            passiveRange = pi/2;
            s.visual.grid_spacing = [-1  0  1]*passiveRange;
            
            %%%
            %%%%%%
            % Define system physics
            s.physics.fluid_density = 1;
            s.physics.passive = 1;
            s.physics.k = 0.074;
            s.physics.b = 0.01;
            s.physics.metabolicRate = 0.05;
            s.physics.maxAcc = 200;
            
            %'speed', 'efficiency', or 'metabolism'
            s.passiveObjectiveFunction = 'speed';
            %s.physics.interaction = 'off';
            s.physics = rmfield(s.physics,'passive')
 
            %Functional Local connection and dissipation metric


            s.A = @(alpha1,alpha2) Inertial_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2]);                        % Joint angles

            s.metric = @(alpha1,alpha2) nLinksHead*Inertial_energy_metric(s.geometry,s.physics,[alpha1,alpha2]);
            %s.metric = @(alpha1,alpha2) eye(2);
          
                    
			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-2*pi,2*pi,-pi,pi];
            
			%densities for various operations
			s.density.vector = [21 21]; %density to display vector field
			s.density.scalar = [51 51]; %density to display scalar functions
			s.density.eval = [31 31];   %density for function evaluations
            s.density.metric_eval = [31 31]; %density for metric evaluation
            s.density.mass_eval = [31 31]; % density for mass matrix evaluation
            s.density.coriolis_eval = [31 31];
            s.density.finite_element=31;

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*pi;
			s.tic_locs.y = [-1 0 1]*pi;


            % System type flag
            s.system_type = 'inertia';
			%%%%
			%Save the system properties
			output = s;


	end

end

