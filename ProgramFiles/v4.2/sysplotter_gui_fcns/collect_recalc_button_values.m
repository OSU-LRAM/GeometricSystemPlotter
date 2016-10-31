function values = collect_recalc_button_values(handles)
% Collect the radio-button values from the recalculation-control interface


	% List the individual blocks in the recalculation-control
	blocks = {'sys_init','shch_init','sys_calc'...
		,'shch_calc','override'};
	
	% Prime the array of values
	values = zeros(size(blocks));
	
	% Loop over the blocks
	for i = 1:length(blocks)
		
		% Get the selected object in that block
		selected_handle = get(handles.([blocks{i} '_group']),'SelectedObject');
		
		% Get the tag associated with that object
		selected_tag = get(selected_handle,'Tag');
		
		% Parse the tag
		values(i) = [1 2 3 4] *...
			strcmp(selected_tag(end-4:end),{'_auto';'force';'block';'idual'}); % Last is "individual"
		
	end



end