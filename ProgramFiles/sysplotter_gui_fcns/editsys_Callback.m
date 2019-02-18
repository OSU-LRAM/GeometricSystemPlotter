% --- Executes on button press in editsys.
function editsys_Callback(hObject, eventdata, handles)
% hObject    handle to editsys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %Get the system and shape change info
    system_index = get(handles.systemmenu,'Value');
    system_names = get(handles.systemmenu,'UserData');
    
    current_system = system_names{system_index};
    
	if strncmp('sysf_',current_system,5)
		
        %Launch the Configuration GUI
        config_func = str2func(current_system); %Get config function for current system
        [matFilePath] = config_func('savepath'); %Get the config savefile path
		sysf_config_gui(current_system,matFilePath);%Call thegui
        

		
	end


end