function set_plot_resolution_from_file(handles)
% Set the plot resolution based on the system file

	% declare the data directory
	configfile = './sysplotter_config';
    pathnames = load(configfile);

	%Get the system and shape change info
    system_index = get(handles.systemmenu,'Value');
    system_names = get(handles.systemmenu,'UserData');
    
    sys = system_names{system_index};

	% Get the setup configuration file and the path to the system files
	configfile = './sysplotter_config';
	load(configfile,'syspath')
	
	% Run the system initialization to get the specified resolutions
	s = absolute_feval(fullfile(syspath, sys),'initialize',pathnames); 
		
	% Set the grid range
	set(handles.vectorresolution,'UserData',s.grid_range);
	set(handles.scalarresolution,'UserData',s.grid_range);
	
	% If the system is something other than 'default', and the "Hold
	% resolution" checkbox is not set,
	% get the vector and scalar resolution and put them in the gui boxes
	if ~any(strcmp(sys,{'default','start_recent','end_recent'}))...
			&& ~get(handles.holdres,'Value')



		
		% Get the string form of the resolutions
		vector_resolution = num2str(s.density.vector);
		scalar_resolution = num2str(s.density.scalar);
		
		% Set the resolution fields
		set(handles.vectorresolution,'String',vector_resolution);
		set(handles.scalarresolution,'String',scalar_resolution);
		
	end
	


	
end