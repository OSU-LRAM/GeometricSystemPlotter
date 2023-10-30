function output = sysf_serpenoid_lowRe_hf_AmpPhase(input_mode,pathnames)
% System file for a low Reynolds

    % Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Viscous swimmer: Serpenoid high frequency AmpPhase'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/ContinuousBackbone/',...
                'Physics/LowReynoldsRFT/'});

		case 'initialize'

            %%%%%%%%%
            % System definition

            %%%
            % Define the geometry of the system

            % this system's shape is defined by a set of basis functions
            % whose weighted sum is curvature as a function of distance
            % along the backbone
            s.geometry.type = 'general curvature';                

            % The specific basis functions are those for a serpenoid curve,
            % in which the curvature varies sinusoidally along the length
            % of the body. This function is normalized such that the body length
            % is taken as one unit long. For this system, we use a
            % wavelength equal to the body length.
            curvdef_name = 'AmpPhase_wave_two_period';
            n_waves = 2;
            AmpPhasecurvature = AmpPhase_wave_generator;
            AmpPhasecurvature = subs(AmpPhasecurvature,'n',n_waves);
            curvdef_parameters = {'A' 'phi'};
            curvdef_fun = make_curvdef(AmpPhasecurvature,curvdef_parameters,curvdef_name,[1,0,0,0],pathnames.syspath,pathnames.sysplotterpath);

            s.geometry.function = curvdef_fun;%{@(s)serpenoid_1(s,n_waves);@(s)serpenoid_2(s,n_waves)};

            % Total length of the swimmer, in real units
            s.geometry.length = 1;
            
            % base the system off of its center frame
            s.geometry.baseframe = 'com-mean';%'com-mean';

            %%%
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1 0 1]*6 %[-1 -0.5 0 0.5 1]*6;
            
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
            
            % Processing details

            %Range over which to evaluate connection
            s.grid_range = [-1,1,-1,1]*12;

            %densities for various operations
            s.density.vector = [11 11]; %density to display vector field
            s.density.scalar = [21 21]; %density to display scalar functions
            s.density.eval = [21 21];   %density for function evaluations
            s.density.metric_eval = [1 1]*11;
            s.density.finite_element=31;

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


function AmpPhaseWave = AmpPhase_wave_generator
% Generate the curvature profile for a triangle wave




    syms A phi omega s n

    % Wave with sin and cosine amplitudes
    AmpPhaseWave = (A*cos(2*pi*s*n - phi));

    % % Amplitude of the wave
    % amp_wave = (a1^2+a2^2)^(.5);
    % 
    % % Normalize the wave
    % norm_wave = raw_wave/amp_wave;
    % 
    % % Pinch the wave
    % pinched_wave = (abs(norm_wave)^n) * sign(norm_wave);
    % 
    % % Restore the scale of the wave
    % pinched_wave_restored = pinched_wave*amp_wave;
    % 
    % % Scale to match max angle with original wave
    % triangle_wave_matched = pinched_wave_restored * (P*[n^3;n^2;n;1]);

end