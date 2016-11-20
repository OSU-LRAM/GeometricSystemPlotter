
% --- Executes on selection change in stretchmenu.
function stretchmenu_Callback(hObject, eventdata, handles,active)
	% hObject    handle to stretchmenu (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Hints: contents = cellstr(get(hObject,'String')) returns stretchmenu contents as cell array
	%        contents{get(hObject,'Value')} returns selected item from stretchmenu
	
	% Update the run display
	refresh_runinfo_Callback([], eventdata, handles);

end
