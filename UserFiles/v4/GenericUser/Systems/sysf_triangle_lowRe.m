function output = sysf_triangle_lowRe(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
    end
    
    
    %%%%
    % Name the file that will hold the curvature definition for this system
    % (its contents are the first elements of the initialization code)
    curvdef_name = 'triangle_wave_one_period';
	
	
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Viscous swimmer: Triangle'; % Display name

		case 'dependency'

            % Dependent on the functions for generating the backbone,
            % connection, and metric from the curvature, and on the curvdef
            % file generated during initialization.
			output.dependency = [fullfile(pathnames.sysplotterpath,...
                {'Utilities/curvature_mode_toolbox/triangle_wave_generator';...
                'Utilities/curvature_mode_toolbox/backbone_from_general_curvature.m';...
				'Utilities/LowRE_toolbox/LowRE_dissipation_metric_from_general_curvature.m';...
				'Utilities/LowRE_toolbox/LowRE_local_connection_from_general_curvature.m'});
                fullfile(pathnames.syspath,curvdef_name,'.mat')];

		case 'initialize'
            
            %%%%
            % Make a curvdef function describing the shape of the system
            
            % Generate a function that uses sin^power to approximate a
            % triangle wave, using the generator function in the curvature
            % tool box (simpler functions could be specified by a string or
            % symbolic expression right here)
            triangle_wave = absolute_feval(fullfile(pathnames.sysplotterpath, 'Utilities/curvature_mode_toolbox/triangle_wave_generator'));
            
            % Substitute in a unit value for the frequency, and 10 for the
            % sharpness power
            triangle_wave = subs(triangle_wave,{'omega','n'},{1 10});
            
            % Use the triangle_wave expression to generate a curvdef
            % file that numerically integrates the curvature and finds
            % the derivative of the curvature with respect to the shape
            % variables. (this file is used by the connection and
            % metric calculations)
            make_curvdef(triangle_wave,{'a1' 'a2'},curvdef_name,[1,0,0,0],pathnames.sysplotterpath,pathnames.syspath)
            
            % Declare this function as the curvdef for the connection and
            % metric calculations
            curvdef_fun = str2func(['curv_' curvdef_name]);
            
            %%%%%%%%%%%
            %%%%%%%%%%%

			%Functional representation of local connection. 
			s.A_num = @(a1,a2) Conn_num(a1,a2,curvdef_fun);
			% power metric
			s.metric = @(a1,a2) M_helper(a1,a2,curvdef_fun);
            
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
			s.density.scalar = [1 1]*11; %density to display scalar functions
			s.density.eval = [1 1]*11;   %density for function evaluations
            s.density.metric_eval = [1 1]*11;
			s.finite_element_density = 11;

			%%%
			%Display parameters

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*6;
			s.tic_locs.y = [-1 0 1]*6;


			%%%%
			%Save the system properties
			output = s;


	end

end

function [A]=Conn_num(a1,a2,curvdef_fun)

	% Get the local connection for the specified shape, with unit length
	A = LowRE_local_connection_from_general_curvature(curvdef_fun,[a1;a2],1,1,2);


end



function M = M_helper(a1,a2,curvdef_fun)

	% Get the local connection for the specified shape, with unit length
	M = LowRE_dissipation_metric_from_general_curvature...
				(curvdef_fun,[a1;a2],1,1,2);

end

