function initialize_filelist_menu(hObject,heading,dirname,prefix,recent_listname,max_recent,handles)
%
    % Load the list of systems and their shortnames
    [displaynames, filenames] = list_files(dirname,prefix);
    
    % Put default and null options on system list, and pad the shortnames
    
    displaylist = [vertcat(heading.display(:));'-';'-';displaynames];
    filelist = [vertcat(heading.file(:));'start_recent';'end_recent';filenames];
        
    % insert displaynames as the string and shortnames as the userdata
    set(hObject,'String',displaylist,'UserData',filelist,'Value',1);
	
	% update for any recent items
	menu_recent_update(hObject,recent_listname,max_recent,1,handles);
	
	% Set the menu to its top entry
	set(hObject,'value',1);
	
end