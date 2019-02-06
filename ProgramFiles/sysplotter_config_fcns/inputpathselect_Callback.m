% --- Executes on button press in systempathselect.
function inputpathselect_Callback(hObject, eventdata, handles)
% hObject    handle to systempathselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Set the directory string to a user-selected path
    oldtargetdir = get(handles.inputpathconfig,'String');
    newtargetdir = uigetdir(get(handles.inputpathconfig,'String'),'Select directory of configuration files');
    
    if ~isnumeric(newtargetdir)
        targetdir = newtargetdir;
    else
        targetdir = oldtargetdir;
    end
    
    set(handles.inputpathconfig,'String',targetdir);
    
    % Verify that the target directory has the necessary subdirectories
    v = verify_configdir(targetdir);
    
    % Only activate the OK button if the string is valid
    if v
        set(handles.OKbutton,'enable','on');
    else
        set(handles.OKbutton,'enable','off');
    end

end