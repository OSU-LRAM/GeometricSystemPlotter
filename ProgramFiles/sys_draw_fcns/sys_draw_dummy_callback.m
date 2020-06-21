function sys_draw_dummy_callback(hObject, eventdata, plot_structure,sys,shch,stretch)
%strip out extra callback info

	handles = guidata(hObject);
	
	%%%%%%%%%%%%%%
	% Get the desired vector and scalar plotting resolutions
	resolution.vector = str2num(get(handles.vectorresolution,'String'));
	resolution.scalar = str2num(get(handles.scalarresolution,'String'));
	resolution.vector_range = get(handles.vectorresolution,'UserData');
	resolution.scalar_range = get(handles.scalarresolution,'UserData');

    sys_draw(plot_structure,sys,shch,stretch,handles.progresspanel,0,resolution,handles);

end