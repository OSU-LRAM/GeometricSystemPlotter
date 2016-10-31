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
			
			%Vectorize the connection definition functions if necessary (should
			%only be necessary if connection is defined as an inline function)
			s = vectorize_connection(s);
			
			%Create grids for evaluating the connection functions
			s = create_grids(s);
			
			%Evaluate the connection and metric over the fine grid for calculations and the coarse
			%grid for vector display
			s = evaluate_connection(s);
			s = evaluate_metric(s);
			
			%Reshape the evaluated connection and metric to format compatible with plotting
			%over base space
			s = reshape_connection(s);
			s = reshape_metric(s);
			
			%Merge components of evaluated connection and metric
			s = merge_connection(s);
			s = merge_metric(s);
			
			%Calculate optimal coordinate choice
			s = optimize_coordinate_choice(s);
			
			%Calculate the height functions from the connection
			s = calc_height_functions(s);
						
			%Save out the updated system properties
			save(outfile,'s')
			
	end
	
end