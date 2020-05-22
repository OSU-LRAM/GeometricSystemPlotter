% --- Executes on button press in any plotpushbutton.
function plot_info = plotpushbutton_Callback_optimize(hObject, eventdata, handles)
% hObject    handle to plotpushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Update the waitbar to indicate that the process has started
waitbar2a(0,handles.progresspanel,'Gathering plot data');

%Remove old subplots
delete(get(handles.plot_thumbnails,'Children'))

%extract name of source
source_name = get(hObject,'Tag');

%extract the column number as text
source_number_text = source_name(end);

%get the checkbox values
% [box_names, box_active, box_values, box_enabled, plot_types, ...
% 	plot_subtypes,merged_plot_subtypes, plot_style] =...
% 	get_box_values(source_number_text,handles); %#ok<ASGLU>

[box_names,...
    box_active,...
    box_values,...
    box_enabled,...
    plot_types, ...
	plot_subtypes,...
    plot_style,...
    plot_coordinates] =...
        get_box_values(source_number_text,handles); %#ok<ASGLU>

%get the CCF type to plot
CCFtype = get(findobj(handles.(['CCFradio' source_number_text]),'Value',1),'Tag');
CCFtype(3) = [];
    
%get the height function type to plot
hfuntype = get(findobj(handles.(['CCFradio' source_number_text]),'Value',1),'Tag');
hfuntype(3) = [];

% Get the state of the Stretch menu (coordinate conversion to flatten
% metric)
stretchstate = get(handles.stretchmenu,'Value');
stretch_names = get(handles.stretchmenu,'UserData');
current_stretch = stretch_names(stretchstate);


	
% Initialize the plot windows
plots_to_make = initialize_plot_windows(box_names,...
                                                 box_active,...
                                                 plot_style,...
                                                 stretchstate,...
                                                 current_stretch,...
                                                 plot_types,...
                                                 plot_subtypes,...
                                                 plot_coordinates,...
                                                 handles,...
                                                 source_number_text,...
                                                 CCFtype);
% plots_to_make = initialize_plot_windows(box_names,box_active,plot_types,merged_plot_subtypes...
% 	,plot_style,hfuntype,stretchstate,handles,source_number_text);


%%%%%%%%%%%
%Get the system and shape change info
system_index = get(handles.systemmenu,'Value');
system_names = get(handles.systemmenu,'UserData');

current_system = system_names{system_index};

shch_index = get(handles.shapechangemenu,'Value');
shch_names = get(handles.shapechangemenu,'UserData');

current_shch = shch_names{shch_index};

%%%%%%%%%%%%%%
% Get the desired vector and scalar plotting resolutions
resolution.vector = str2num(get(handles.vectorresolution,'String'));
resolution.scalar = str2num(get(handles.scalarresolution,'String'));
resolution.vector_range = get(handles.vectorresolution,'UserData');
resolution.scalar_range = get(handles.scalarresolution,'UserData');


% current_shch='null';
%Call the draw function
%    test_plot(plots_to_make,current_system,current_shch)
plot_info = sys_draw(plots_to_make,current_system,current_shch,current_stretch,handles.progresspanel,1,resolution,handles);

%Show full progress bar
waitbar2a(0.8,handles.progresspanel,'Optimizing')


end