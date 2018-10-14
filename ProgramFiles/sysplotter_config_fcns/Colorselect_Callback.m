% --- Executes on button press in Refpointselect.
function Colorselect_Callback(hObject, eventdata, handles)
% hObject    handle to Refpointselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Set the directory string to a user-selected path
    [targetfile, targetpath, targetindex] = uigetfile('color_*.m','Select file with color specification', get(handles.Colorconfig,'String'));
    
    try
    addpath(targetpath);
    end
    
    
    set(handles.Colorconfig,'String',fullfile(targetpath,targetfile));
	

end
