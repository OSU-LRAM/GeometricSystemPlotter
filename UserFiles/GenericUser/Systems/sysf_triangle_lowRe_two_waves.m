function output = sysf_triangle_lowRe_two_waves(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
    end
    
    if ~exist('pathnames','var')
        
        pathnames = load('sysplotter_config');
        
    end
    
    
    %%%%
    % Name the file that will hold the curvature definition for this system
    % (its contents are the first elements of the initialization code)
    curvdef_name = 'triangle_wave_two_periods';
	
	
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Viscous swimmer: Triangle two waves'; % Display name

		case 'dependency'

            % Dependent on the functions for generating the backbone,
            % connection, and metric from the curvature
            output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/ContinuousBackbone/',...
                'Physics/LowReynoldsRFT/'});


		case 'initialize'
            
            %%%%
            % Make a curvdef function describing the shape of the system
            
            % Generate a function that uses sin^power to approximate a
            % triangle wave, using the generator function below (simpler
            % functions could be specified by a string or 
            % symbolic expression right here)
            triangle_wave_curvature = triangle_wave_generator;
            
            % Substitute in a unit value for the frequency, and 10 for the
            % sharpness power
            triangle_wave_curvature = subs(triangle_wave_curvature,{'omega','n'},{2 20});
            
            % Use the triangle_wave expression to generate a curvdef
            % file that numerically integrates the curvature and finds
            % the derivative of the curvature with respect to the shape
            % variables. (this file is used by the connection and
            % metric calculations)
            curvdef_parameters = {'a1' 'a2'};
            curvdef_fun = make_curvdef(triangle_wave_curvature,curvdef_parameters,curvdef_name,[1,0,0,0],pathnames.syspath,pathnames.sysplotterpath);
             
            %%%%%%%%%%%
            %%%%%%%%%%%

            %%%%%%%%%
            % System definition

            %%%
            % Define the geometry of the system

            % this system's shape is defined by the general curvature
            % function created above
            s.geometry.type = 'general curvature';                
            s.geometry.function = curvdef_fun;
            s.geometry.baseframe = 'center';
            
            % Total length of the swimmer, in real units
            s.geometry.length = 1;
            
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = ndgrid([-1 -0.5 0 0.5 1]*6+.001); %.001 is to avoid a singularity at zero

            
            %%%

            %%%
            % Define the physics of the system

            % This system is treated as moving in a viscous fluid. Key elements we
            % need for this model are a drag coefficient and a drag ratio
            s.physics.drag_coefficient = 1;                   % multiplier from longitudinal velocity to drag force. Changing this scales the dissipation matrix
            s.physics.drag_ratio = 2;                         % ratio of lateral:longitudinal drag

            
             % Locomotion model derived from viscous drag forces reacting to
            % local velocities of elements on the body
            s.A = @(a1,a2) LowRE_local_connection(s.geometry,s.physics,[a1;a2]);
            
            % Locomotion model derived from viscous drag forces reacting to
            % local velocities of elements on the body
            s.metric = @(a1,a2) LowRE_dissipation_metric(s.geometry,s.physics,[a1;a2]);
            %%%

                        
			%%%
			%Processing details

            % There is a singularity at the origin (because the raw wave is
            % divided by its amplitude when making it triangular). This
            % field suppresses a warning message about the singularity.
            s.ignore_singularity_warning = 1;
            
            
			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*12;
            
			%densities for various operations
			s.density.vector = [1 1]*11; %density to display vector field
			s.density.scalar = [1 1]*51; %density to display scalar functions
			s.density.eval = [1 1]*31;   %density for function evaluations
            s.density.metric_eval = [1 1]*11;
            s.density.finite_element=31;

			%%%
			%Display parameters

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*6;
			s.tic_locs.y = [-1 0 1]*6;

            % Set system type variable for gait optimization
            s.system_type = 'drag';
            
			%%%%
			%Save the system properties
			output = s;


	end

end

% function [A]=Conn_num(a1,a2,curvdef_fun)
% 
% 	% Get the local connection for the specified shape, with unit length
% 	A = LowRE_local_connection_from_general_curvature(curvdef_fun,[a1;a2],1,1,2);
% 
% 
% end
% 
% 
% 
% function M = M_helper(a1,a2,curvdef_fun)
% 
% 	% Get the local connection for the specified shape, with unit length
% 	M = LowRE_dissipation_metric_from_general_curvature...
% 				(curvdef_fun,[a1;a2],1,1,2);
% 
% end

function triangle_wave_matched = triangle_wave_generator
% Generate the curvature profile for a triangle wave

    % Generate a scale factor so that the maximum angle reached by the sin and
    % pinched sin basis functions is the same
    n = logspace(-1,1.5,100);

    x = linspace(0,pi/2);

    base_area = trapz(x,cos(x));

    for idx = 1:numel(n)

        scale_factor(idx) = base_area/trapz(x,cos(x).^n(idx)); %#ok<AGROW>

    end

    P = polyfit(n,scale_factor,3);



    % Now generate the pinched sin function


    syms a1 a2 omega s n

    % Wave with sin and cosine amplitudes
    raw_wave = (a1*cos(omega*2*pi*s)+a2*sin(omega*2*pi*s));

    % Amplitude of the wave
    amp_wave = (a1^2+a2^2)^(.5);

    % Normalize the wave
    norm_wave = raw_wave/amp_wave;

    % Pinch the wave
    pinched_wave = (abs(norm_wave)^n) * sign(norm_wave);

    % Restore the scale of the wave
    pinched_wave_restored = pinched_wave*amp_wave;

    % Scale to match max angle with original wave
    triangle_wave_matched = pinched_wave_restored * (P*[n^3;n^2;n;1]);

end


