% --- Executes on button press in gait_gui_draw.
function gait_gui_draw_Callback(hObject, eventdata, handles)
% hObject    handle to gait_gui_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the last plot pushbutton used
if isfield(handles.figure1.UserData,'lastpushbutton')
    lastpushbutton = handles.figure1.UserData.lastpushbutton;
else
    lastpushbutton = 'plotpushbutton1';
end

% Plot whatever gaits are selected in the third checkbox column
plot_info = plotpushbutton_Callback(findall(0,'tag',lastpushbutton), eventdata, handles);

% Execute the gait_gui_draw command
    gait_gui_draw(plot_info(1).axes(1),hObject, eventdata, handles);
    
%     set(handles.shapechangemenu,'Value',rn(1));
%     plot_info = plotpushbutton_Callback(hObject, eventdata, handles); 

end
