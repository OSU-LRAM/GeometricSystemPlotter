function guiDisplaySelection_Callback(hObject, eventdata, handles)
% hObject    handle to guiDisplayOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns guiDisplayOptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from guiDisplayOptions

propertyMenuValue = get(handles.guiDisplaySelection,'Value');
% % GUIposition = get(handles.figPendulum,'Position');
%
if propertyMenuValue == 1 % default
    if ispc
        propertyFile = 'propertyDataWindows.mat';
    elseif ismac
        propertyFile = 'propertyDataMacOS.mat';
    elseif isunix
        propertyFile = 'propertyDataLinux.mat';
    end
elseif propertyMenuValue == 2 % PC
    propertyFile = 'propertyDataWindows.mat';
elseif propertyMenuValue == 3 % MAC
    propertyFile = 'propertyDataMacOS.mat';
elseif propertyMenuValue == 4 % MAC
    propertyFile = 'propertyDataLinux.mat';
elseif propertyMenuValue == 5 % Custom (currently unavailable)
    %         uiopen
    propertyFile = uigetfile('*.mat','Select Property File');
end

set(handles.displayConfigFile,'String',propertyFile)

end