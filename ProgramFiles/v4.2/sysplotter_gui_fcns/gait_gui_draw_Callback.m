% --- Executes on button press in gait_gui_draw.
function gait_gui_draw_Callback(hObject, eventdata, handles)
% hObject    handle to gait_gui_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot whatever gaits are selected in the third checkbox column
plot_info = plotpushbutton_Callback(findall(0,'tag','plotpushbutton3'), eventdata, handles);

% Execute the gait_gui_draw command
    gait_gui_draw(plot_info(1).axes(1));

end
