% --- Executes on button press in OKbutton.
function OKbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OKbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Save the system path string to the configuration data file
    inputpath = get(handles.inputpathconfig,'String'); %#ok<NASGU>
	% path to generated data

	% Get the version number of sysplotter and use it to make a datapath
	currentDirectory = pwd;
	[~, deepestFolder, ~] = fileparts2(currentDirectory);
	syspath = fullfile(inputpath,'Systems');
	shchpath = fullfile(inputpath,'Shape_Changes');
	stretchpath = fullfile(inputpath,'Stretches');
	datapath = fullfile(inputpath,  '/sysplotter_data/');
    
    % Get the colors to use in plots
    [~,colorfunction, ~ ] = fileparts2(get(handles.Colorconfig,'String'));
    Colorset = feval(colorfunction);
    Colorpath = get(handles.Colorconfig,'String');
    
    % Additional hard-coded paths
    sysplotterpath = pwd;
    
    % propertyfile
    propertyfilepath = fullfile(sysplotterpath,get(handles.displayConfigFile,'String'));

    % Save the path info to a file for sysplotter to refer to
    save('sysplotter_config','inputpath','syspath','shchpath','stretchpath','datapath','Colorset','Colorpath','sysplotterpath','propertyfilepath');
    
    % Update the sysplotter_inputpath variable in the workspace
    assignin('base','sysplotter_inputpath',inputpath);
        

    
    % Close the window
    delete(handles.sysplotter_config_gui);

end