
% --- Executes when selected object is changed in override_group.
function override_group_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in override_group 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

	% List the individual blocks in the recalculation-control
	blocks = {'sys_init_group','shch_init_group','sys_calc_group','shch_calc_group'};
	
	% turn block names into handles
	block_handles = cellfun(@(b) handles.(b),blocks,'UniformOutput',false);
    block_handles = [block_handles{:}];
	
	% block children
	block_children = [block_handles(:).Children];
	
	% block children that are radio handles
	block_radio = block_children(strcmp('radiobutton',get(block_children,'Style')));

	% If the "Individual" button is selected, activate all of the section
	% blocks. Otherwise, disable the blocks
	if eventdata.NewValue == handles.override_individual
		
		set(block_radio,'Enable','on')
		
	else
		
		set(block_radio,'Enable','off')
		
	end
		
		
	% Update the run display
	refresh_gui_Callback([], eventdata, handles);


end