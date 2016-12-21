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
			
			%Create grids for evaluating the connection functions
			s = create_grids(s);
			
            % Ensure that there is a metric
            s = ensure_connection_and_metric(s);

            
            
			%Evaluate the connection and metric over the fine grid for calculations and the coarse
			%grid for vector display
			s = evaluate_connection(s);
            s = evaluate_metric(s);
			
			
			%Merge components of evaluated connection and metric
			s = merge_connection(s);
			s = merge_metric(s);
			
			%Calculate optimal coordinate choice
			s = optimize_coordinate_choice(s);
			
			%Calculate the height functions from the connection
			s = calc_constraint_curvature(s);
            
            %Build a stretch function corresponding to the metric
            s = calc_stretch_functions(s);
						
			%Save out the updated system properties
			save(outfile,'s')
			
	end
	
end