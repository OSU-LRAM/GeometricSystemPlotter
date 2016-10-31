function resolution_box_Callback(hObject, eventdata, handles)
% hObject    handle to calling resolution box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vectorresolution as text
%        str2double(get(hObject,'String')) returns contents of vectorresolution as a double


% If the text is hand-edited in the resolution box, set the checkbox to
% hold the specified resolution even if the system is changed
set(handles.holdres,'Value',1);

end