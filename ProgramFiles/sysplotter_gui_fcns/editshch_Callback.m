% --- Executes on button press in editshch.
function editshch_Callback(hObject, eventdata, handles)
% hObject    handle to editshch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    shch_index = get(handles.shapechangemenu,'Value');
    shch_names = get(handles.shapechangemenu,'UserData');
    
    current_shch = shch_names{shch_index};
	
	if strncmp('shchf_',current_shch,6)
		
		edit([current_shch '.m']);
		
	end

end