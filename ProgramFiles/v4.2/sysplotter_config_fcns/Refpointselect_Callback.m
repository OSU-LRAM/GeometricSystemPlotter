% --- Executes on button press in Refpointselect.
function Refpointselect_Callback(hObject, eventdata, handles)
% hObject    handle to Refpointselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Set the directory string to a user-selected path
    targetdir = uigetdir(get(handles.Refpointconfig,'String'),'Select directory of configuration files');
    set(handles.Refpointconfig,'String',targetdir);
	

end
