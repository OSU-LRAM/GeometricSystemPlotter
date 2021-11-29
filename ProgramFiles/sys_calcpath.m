function [output] = sys_calcpath(input_mode,systemfilename,shchfilename)
%Calculate the effects of running the path defined in pathfilename on the
%system defined in systemfilename

	%filename for system data with path information
	systemshchfilename = [systemfilename '__' shchfilename];

	% Get the setup configuration file
	configfile = './sysplotter_config';
	load(configfile,'datapath')
	
	infile1 = fullfile(datapath, [systemfilename '_calc.mat']);
	infile2 = fullfile(datapath, [shchfilename '.mat']);
	outfile = fullfile(datapath, [systemshchfilename  '.mat']);

	switch input_mode
		
		case 'dependency'
			
			output.dependency = {'sys_calcpath_fcns/';
				infile1;
				infile2};
			output.product = outfile;
			
		case 'calculate'	
						
			output = systemshchfilename;
			
			%Load the system properties from the system data file
			load(infile1,'s')
			
			%convert empty matrix input to 'null' string for shchfilename
			if isempty(shchfilename)
				
				shchfilename = 'null';
				
			end
			
			
			%If there's a path to evaluate, do so
			if ~strcmp('null',shchfilename)
				
				%load the path properties from the path data file
				load(infile2,'p')
				
				%Vectorize the shape change equation if it's an inline function
				p = vectorize_shch(p); %#ok<NODEF>
                
				
				%get shape change loci
				p = find_loci(s,p);
				
				%Get the resulting fiber motion
				p = find_g_sys(s,p); %#ok<NASGU>
				
				
				
				%Otherwise, save out a data file for plotting a system with no path
				%included
			else
				
				%file
				p = struct([]); %#ok<NASGU>
				
			end
						
			save(outfile,'s','p')
			
	end
    
end


    


