% --- Executes when selected object is changed in sys_init_group.
function sys_init_group_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in sys_init_group 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% Update the run display
	refresh_gui_Callback([], eventdata, handles);

end