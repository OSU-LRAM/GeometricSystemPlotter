% --- Executes during object creation, after setting all properties.
function stretchmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stretchmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    set(hObject,'BackgroundColor','white');
    
    
    % Set the menu options to be unstretched or metric-stretched
    set(hObject,'String',{'No stretch';'Metric stretch';'Metric surface'});
    

  
end