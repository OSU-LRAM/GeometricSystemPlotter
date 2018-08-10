% --- Executes on button press in animatebutton_gui2.
function animatebutton_gui2_Callback(hObject, eventdata, handles)
% hObject    handle to animatebutton_gui2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the number of gaits
info_needed.Number_gaits = str2double(get(handles.number_input,'string'));

% Get the Frame rate
info_needed.Framerate = str2double(get(handles.framrate,'string'));



%Get the system info
system_index = get(handles.old.systemmenu,'Value');
system_names = get(handles.old.systemmenu,'UserData');

current_system = system_names{system_index};

% Remove the "sysf_" from the system
current_system2 = current_system(6:end);

%%%%%%%%%%%%%%%% Select an animate system function

%Load system to check if there's a BackboneShape field
configfile = './sysplotter_config';
load(configfile,'datapath')
sysfile = fullfile(datapath, current_system);
load(sysfile,'s')
    
% Check whether it is 3-links, 4-inks, serpenoid ot triangular
if findstr(current_system2,'serpenoid')
    
    animate_system = @animate_serpenoid_swimmer;
    
elseif findstr(current_system2,'triangular')
    
    animate_system = @animate_triangular_swimmer;
    
elseif findstr(current_system2,'four_link')
    
    animate_system = @animate_4_links_swimmer;
    
elseif findstr(current_system2,'floating')
    
    animate_system = @animate_floating_snake;
    
elseif isfield(s,'BackboneShape')
% If there isn't a specialized animation function but the system has a
% geometry specified, use animate_backbone

    animate_system = @animate_backbone;

else
    
    animate_system = @animate_3_links_swimmer;
    
end

%%%%%%%%%%%%%%%%

% Get the shape change info
shch_index = get(handles.old.shapechangemenu,'Value');
shch_names = get(handles.old.shapechangemenu,'UserData');

current_shch = shch_names{shch_index};

% Remove "shchf_" from thr shape change
current_shch2 = current_shch(7:end);

% Find the path to user file
current_directory = pwd;
info_needed.current_directory = current_directory;

info_needed.UserFile_path = datapath;

% Passed to animation() for use in movie export filename
info_needed.current_system2 = current_system2;
info_needed.current_shch2 = current_shch2;

% Run the animation
animate_system([],info_needed)
end
