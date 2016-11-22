% --- Executes on button press in HodgeHelmholtzselect.
function HodgeHelmholtzselect_Callback(hObject, eventdata, handles)
% hObject    handle to HodgeHelmholtzselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Set the directory string to a user-selected path
    targetdir = uigetdir(get(handles.HodgeHelmholtzconfig,'String'),'Select directory of configuration files');
    set(handles.HodgeHelmholtzconfig,'String',targetdir);



end
