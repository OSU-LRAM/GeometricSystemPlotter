% --- Executes on selection change in systemmenu.
function systemmenu_Callback(hObject, eventdata, handles)
    % hObject    handle to systemmenu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns systemmenu contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from systemmenu
    
    
	% Maximum number of recent items to keep (note that max_recent is set
	% in systemmenu_CreateFcn when sysplotter is initialized -- to most
	% effectively change the memory depth, edit both values).
	max_recent = 5;
	
	% Update the menu recent items
	menu_recent_update(hObject,'system_list',max_recent,1,handles);
	
	% Trigger an update of the shape change list
	shapechangemenu_Callback(handles.shapechangemenu,eventdata,handles,0);
		
	% Trigger an update of the shape change list
	stretchmenu_Callback(handles.stretchmenu,eventdata,handles,0);

	%Disable the plotting button unless a reasonable value has been set in
    %this menu
    plotpublishpushbutton_enable_menu(handles)
	
	% If a system was selected, set the resolution fields to the natural
	% resolution for the system, unless the field is locked
    set_plot_resolution_from_file(handles);


end