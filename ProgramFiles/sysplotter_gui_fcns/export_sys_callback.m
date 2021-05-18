% --- Executes on button press in pushbutton22
function [s,p] = export_sys_callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	% Get the setup configuration file
	configfile = './sysplotter_config';
	load(configfile,'datapath')
    
    %Get system name from dropdown
    system_index = get(handles.systemmenu,'Value');
    system_names = get(handles.systemmenu,'UserData');
    sys = system_names{system_index};
    
    %Get path name from dropdown
    shch_index = get(handles.shapechangemenu,'Value');
    shch_names = get(handles.shapechangemenu,'UserData');
    shch = shch_names{shch_index};
    
    %Get stretch
    stretchstate = get(handles.stretchmenu,'Value');
    stretch_names = get(handles.stretchmenu,'UserData');
    stretch = lower(stretch_names(stretchstate));
    
    %Update system if out of date
    sys_update(sys,shch,stretch,handles.progresspanel,handles)
    
    %Update waitbar
    waitbar2a(.8,handles.progresspanel,'Loading Data');
    
    %merge the system and shape change file names into one
	plot_data_file = [sys '__' shch];
    
    if strcmpi(sys,'default') || strcmpi(shch,'null')
        error('Must specify system and path to export');
    end
    
    %load the system and path data
	load(fullfile(datapath, plot_data_file),'s','p');
    assignin('base','s',s);
    assignin('base','p',p);
    
    waitbar2a(1,handles.progresspanel,'Loaded Data');
end