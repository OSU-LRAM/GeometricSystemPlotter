% --- Executes on button press in refresh_gui.
function refresh_gui_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Refresh the menus and the runinfo gui indicators
    refresh_menus_Callback(hObject, eventdata, handles);
    refresh_runinfo_Callback(hObject, eventdata, handles);



end