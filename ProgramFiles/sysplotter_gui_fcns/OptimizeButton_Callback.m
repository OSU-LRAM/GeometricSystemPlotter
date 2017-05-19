% --- Executes on button press in pushbutton20.
function OptimizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 plot_info = plotpushbutton_Callback_optimize(findall(0,'tag','plotpushbutton3'), eventdata, handles);

% Execute the gait_gui_draw command
    gait_gui_optimize(plot_info(1).axes(1),hObject, eventdata, handles);
waitbar2a(1,handles.progresspanel,'Finished Plotting')
end