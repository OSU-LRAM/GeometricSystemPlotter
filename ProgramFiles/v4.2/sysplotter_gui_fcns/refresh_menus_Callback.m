% --- Executes on button press in refresh_menus.
function refresh_menus_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_menus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% Refreshes the menus

	% gather data from the system file menu
	systemdata = get(handles.systemmenu,'UserData');
	systemvalue = get(handles.systemmenu,'Value');
	
	old_system_setting = systemdata{systemvalue};
	
	% gather data from the shape change menu
	shchdata = get(handles.shapechangemenu,'UserData');
	shchvalue = get(handles.shapechangemenu,'Value');
	
	old_shapechange_setting = shchdata{shchvalue};

	% gather data from the stretch menu
	stretchdata = get(handles.stretchmenu,'UserData');
	stretchvalue = get(handles.stretchmenu,'Value');
	
	old_stretch_setting = stretchdata{stretchvalue};
	
	%%%%%%%%%%%%
	% Update the system menu
	systemmenu_CreateFcn(handles.systemmenu,eventdata,handles);
	
	% Get the new file listing
	systemdata = get(handles.systemmenu,'UserData');
	
	% Find the previous system's location in the new list
	new_system_value = find(strcmp(old_system_setting,systemdata),1);
	if isempty(new_system_value)
		new_system_value = 1;
	end
	set(handles.systemmenu,'value',new_system_value);
	

	%%%%%%%%%%%%
	% Update the shape change menu
	shapechangemenu_CreateFcn(handles.shapechangemenu,eventdata,handles)
	
	% Get the new file listing
	shchdata = get(handles.shapechangemenu,'UserData');
	
	
	% Find the previous shape change's location in the new list
	new_shapechange_value = find(strcmp(old_shapechange_setting,shchdata),1);
	if isempty(new_shapechange_value)
		new_shapechange_value = 1;
	end
	set(handles.shapechangemenu,'value',new_shapechange_value);

	
	%%%%%%%%%%%%
	% Update the stretch menu
	stretchmenu_CreateFcn(handles.stretchmenu,eventdata,handles)
	
	% Get the new file listing
	stretchdata = get(handles.stretchmenu,'UserData');
	
	
	% Find the previous shape change's location in the new list
	new_stretch_value = find(strcmp(old_stretch_setting,stretchdata),1);
	if isempty(new_stretch_value)
		new_stretch_value = 1;
	end
	set(handles.stretchmenu,'value',new_stretch_value);

	
end