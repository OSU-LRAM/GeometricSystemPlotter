
% --- Executes on selection change in stretchmenu.
function stretchmenu_Callback(hObject, eventdata, handles,active)
	% hObject    handle to stretchmenu (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Hints: contents = cellstr(get(hObject,'String')) returns stretchmenu contents as cell array
	%        contents{get(hObject,'Value')} returns selected item from stretchmenu

	% If called externally with active=0, will pass on the information that
	% this call was triggered by something other than menu selection to the
	% recent-updater
	if ~exist('active','var')
		active = 1;
	end
	
	% Maximum number of recent items to keep (note that max_recent is set
	% in shapechangemenu_CreateFcn when sysplotter is initialized -- to most
	% effectively change the memory depth, edit both values).
	max_recent = 5;


	%%%%%%%%%
	% Update the menu recent items
	
	% get the filename of the currently selected system
	all_systems = get(handles.systemmenu,'UserData');
	selected_system = all_systems{get(handles.systemmenu,'Value')};
	
	% Only set the recent items if this is an actual system
	non_system_items = {'default','start_recent','end_recent'};
	if ~any(strcmp(selected_system,non_system_items))
		recent_list_name = [selected_system '_stretch'];
	else
		recent_list_name = [];
	end
	
	% look for any stretches recently used with this system
	menu_recent_update(hObject,recent_list_name,max_recent,active,handles);
	
	% Update the run display
	refresh_runinfo_Callback([], eventdata, handles);

end
