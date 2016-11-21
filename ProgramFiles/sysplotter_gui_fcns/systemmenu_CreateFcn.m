% --- Executes during object creation, after setting all properties.
function systemmenu_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to systemmenu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    %if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    %end
	
	% Load the path files to the system data
	configfile = './sysplotter_config.mat';
	load(configfile);
	
	
    % Maximum number of recent items (note that max_recent is reset
	% in systemmenu_Callback recent items are updated -- to most
	% effectively change the memory depth, edit both values).
    max_recent = 5; 
	
	% Menu heading
	heading.display = {'Select System'};
	heading.file = {'default'};

	% Initialize the menu
	initialize_filelist_menu(hObject,heading,syspath,'sysf','system_list',max_recent,handles);
        
end