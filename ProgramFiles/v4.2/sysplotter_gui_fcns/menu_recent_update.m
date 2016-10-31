function menu_recent_update(hObject,listname,max_recent,active_selection,handles) %#ok<INUSD>
% Update the recent-item list at the top of the menu whose handle is
% hObject. active_selection indicates that the recent_update was being
% triggered by something other than a direct menu selection, and that only
% the contents should be updated, without acting as a selection


	%%%%%%%%%%
	% Place the selected system at the top of the list of recent items (but
	% leave the menu at the actual selected location if on the main list
	
	% Lists of names for objects in menu
	displaylist = get(hObject,'String');
	filelist = get(hObject,'UserData');
	
	% Get the list of filenames already in the recent section of the menu
	start_recent = strmatch('start_recent',filelist);
	end_recent = strmatch('end_recent',filelist);
	extant_recent_filenames = filelist(start_recent+1:end_recent-1);
	
	current_value = get(hObject,'Value');
	current_data = filelist{current_value};
	
	% Load the recent systems file if it exists and a listname is specified
	configfile = './sysplotter_config.mat';
	load(configfile);
	recentfile = [datapath '/list_files_recent.mat'];
	if exist(recentfile,'file') && ~isempty(listname)
		
		load(recentfile); % Load a file listing the most recent systems selected
		
		% If there is already recent-data for this list, load it;
		% otherwise, create an empty list pair
		if isfield(recent_data,listname) %#ok<NODEF>
			recent_displaynames = recent_data.(listname).displaynames;
			recent_filenames = recent_data.(listname).filenames;
			
			% Eliminate any recent names that have since disappeared
			extant_file = cellfun(@(x) any(strcmp(x,filelist)),recent_filenames);
			recent_displaynames = recent_displaynames(extant_file);
			recent_filenames = recent_filenames(extant_file);
			
			
		else
			recent_displaynames = {};
			recent_filenames = {};
		end
		
	else
		
		recent_displaynames = {};
		recent_filenames = {};
		
	end
	
	% Check to make sure this is an actual selection
	non_file_items = {'default','null','start_recent','end_recent'};
	
	if ~isempty(listname) % If there is a list to care about
		
		if ~any(strcmp(filelist{get(hObject,'Value')},non_file_items)) % If the selected item is an actual option
			
			if active_selection % If the menu was selected directly, and is 
				%not updating in response to another

				% Add the selected file to the top of the list of recent files
				new_recent_filenames_redundant = [filelist{get(hObject,'Value')}; recent_filenames];
				new_recent_displaynames = [displaylist{get(hObject,'Value')}; recent_displaynames];

				% Remove any earlier instances of this file in the recent list
				[junk, uniqueI] = unique(new_recent_filenames_redundant,'first'); %#ok<ASGLU>
				new_recent_filenames = new_recent_filenames_redundant(sort(uniqueI));
				new_recent_displaynames = new_recent_displaynames(sort(uniqueI));

				% Trim list to its maximum size
				if length(new_recent_filenames) > max_recent
					new_recent_filenames = new_recent_filenames(1:max_recent);
					new_recent_displaynames = new_recent_displaynames(1:max_recent);
				end
				
			else % Update is in response to an external event, so recent 
				%names should just be pulled from file and current selection should not be added
				new_recent_filenames = recent_filenames;
				new_recent_displaynames = recent_displaynames;
			end
			
		else % selected item is not actually a file (is a header or divider)
			% Just keep the recent names
			
			new_recent_filenames = recent_filenames;
			new_recent_displaynames = recent_displaynames;

		end

	else % There is no list, empty the recent list in the menu
		new_recent_filenames = {};
		new_recent_displaynames = {};
	end
	
	% Check if the list-length has changed size
	lengthdiff = length(new_recent_filenames)-length(extant_recent_filenames);
	
	% Insert the updated recent list into the menu
	
	new_filelist = [filelist(1:start_recent);...
		new_recent_filenames;...
		filelist(end_recent:end)];
	
	new_displaylist = [displaylist(1:start_recent);...
		new_recent_displaynames;...
		displaylist(end_recent:end)];
	
	set(hObject,'String',new_displaylist,'UserData',new_filelist);
	
	% Make sure the selected item stays consistent, even if the number of
	% recent items changes, or the selected item was in the recent items
	% but not at the top
	
	
	% If the listname is empty and this is not an active selection, set the menu to its first entry
	if isempty(listname) && ~active_selection
		
		new_value = 1;
	
	elseif ~active_selection % If the update was not an active selection, set the menu
		% the first appearance of the selected filename in the list
		
		new_value = find(strcmp(current_data,new_filelist),1);
		
	% If the current value was before the recent-list
	elseif current_value <= start_recent
		
		new_value = current_value;
		
	% If the selected item was already in the recent-items list
	elseif (current_value > start_recent) && (current_value < end_recent)
		
		new_value = start_recent + 1;
		
	else %current_value >= end_recent % if current object is after the recent item list
		
		new_value = current_value + lengthdiff;
		
	end
		
	set(hObject,'Value',new_value);
	
	% Save the new recent items out to their data file
	if ~isempty(listname)
		recent_data.(listname).filenames = new_recent_filenames;
		recent_data.(listname).displaynames = new_recent_displaynames;
		save(recentfile,'recent_filenames','recent_data')
	end
	
end