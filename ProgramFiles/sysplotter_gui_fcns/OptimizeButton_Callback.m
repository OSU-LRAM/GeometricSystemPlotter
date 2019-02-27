% --- Executes on button press in pushbutton20.
function OptimizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the last plot pushbutton used
if isfield(handles.figure1.UserData,'lastpushbutton')
    lastpushbutton = handles.figure1.UserData.lastpushbutton;
else
    lastpushbutton = 'plotpushbutton1';
end

    plot_info = plotpushbutton_Callback_optimize(findall(0,'tag',lastpushbutton), eventdata, handles);

% Execute the gait_gui_draw command
    gait_gui_optimize(plot_info(1).axes(1),hObject, eventdata, handles);
    waitbar2a(1,handles.progresspanel,'Finished Plotting')
end