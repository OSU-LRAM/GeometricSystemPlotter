function output = sys_calcsystem(input_mode,systemfilename)
%Numerically evaluate the local connection and height function from the
%symbolic defination, including special accounting for singularities.

	% Get the setup configuration file
	configfile = './sysplotter_config';
	load(configfile,'datapath')
	
	infile = fullfile(datapath, [systemfilename '.mat']);
	outfile = fullfile(datapath, [systemfilename  '_calc.mat']);

	switch input_mode
		
		case 'dependency'
						
			output.dependency = {'sys_calcsystem_fcns';
								infile};
			output.product = {outfile};
			
		case 'calculate'
			
			output = systemfilename;
			
			%Load the system properties from the data file
			load(infile,'s')
			
			%%%%%
			%Processing stages, passed off to secondary functions
            
            % Ensure that there is a connection and a metric
            s = ensure_connection_and_metric(s);
                        
			%Create grids for evaluating the connection functions
			s = create_grids(s);
            
			%Evaluate the connection and metric over the fine grid for calculations and the coarse
			%grid for vector display
			s = evaluate_connection(s);
            s = evaluate_metric(s);
            s = evaluate_metric_derivatives(s);
			
			%Merge components of evaluated connection and metric
			s = merge_connection(s);
			s = merge_metric(s);
			
			%Calculate optimal coordinate choice, respecting group
            optimize = s.conf_space.optimizer_fn;
            s = optimize(s);
			
			%Calculate the constraint curvature functions from the connection
			s = calc_ccf_sys(s);
            
            %Build a stretch function corresponding to the metric only if
            %the system is 2 dimensional
            if length(s.grid_range)/2<3
                s = calc_stretch_functions(s);
            end
            
            % Check for whether system type variable exists in structure; default to
            % 'drag' type if non-existent and throw an error if an unsuitable value is
            % present
            if isfield(s,'system_type')
                if ~( strcmpi(s.system_type,'drag') || strcmpi(s.system_type,'inertia') )
                    error('Invalid system_type flag in system struct; must be ''drag'' or ''inertia''')
                end
            else
                disp('No system_type flag set in system struct, defaulting to ''drag''')
                s.system_type = 'drag';
            end
						
			%Save out the updated system properties
			save(outfile,'s')
			
	end
	
end