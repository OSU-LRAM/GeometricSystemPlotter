% --- Executes during object creation, after setting all properties.
function shapechangemenu_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to shapechangemenu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    set(hObject,'BackgroundColor','white');
    
	% Load the path files to the system data
	configfile = './sysplotter_config.mat';
	load(configfile);

	% Maximum number of recent items (note that max_recent is reset
	% in systemmenu_Callback recent items are updated -- to most
	% effectively change the memory depth, edit both values).
    max_recent = 5; 
	
	% Menu heading
	heading.display = {'No Shape Change'};
	heading.file = {'null'};
	
	% Check the current system menu setting
	if ~isempty(handles)
		system_data = get(handles.systemmenu,'UserData');
		system_value = get(handles.systemmenu,'Value');

		if ~any(strcmp(system_data{system_value},{'default','start_recent','end_recent'}))
			recent_listname = system_data{system_value};
		else
			recent_listname = [];
		end
	else
		recent_listname = [];
	end
	
	% Initialize the menu
	initialize_filelist_menu(hObject,heading,shchpath,'shchf',recent_listname,max_recent,handles);
  
end