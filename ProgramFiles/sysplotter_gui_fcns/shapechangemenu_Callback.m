
% --- Executes on selection change in shapechangemenu.
function shapechangemenu_Callback(hObject, eventdata, handles,active)
    % hObject    handle to shapechangemenu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns shapechangemenu contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from shapechangemenu
    

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
		recent_list_name = selected_system;
	else
		recent_list_name = [];
	end
	
	% look for any shape changes recently used with this system
	menu_recent_update(hObject,recent_list_name,max_recent,active,handles);

	
	
	%Disable the plotting button unless a reasonable value has been set in
    %this menu
    plotpublishpushbutton_enable_menu(handles)

    %enable or disable shape change dependent plots, based on whether or
    %not "No Shape Change" has been selected
    enable_disable_shch_plots(hObject,eventdata,handles)
	
	% Update the run display
	refresh_runinfo_Callback([], eventdata, handles);

end
