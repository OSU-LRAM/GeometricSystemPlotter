% Callback function for the stretch paths. If the string in the box is a
% path to a .mat file containing a convert structure, enable the stretch
% button. If it is not a valid mat file, disable the stretch button, and if
% it is non-empty, color text red to show non-convert structure, or
% background red to show invalid path
function stretchpath_Callback(hObject, eventdata, handles)

	% Add the path to the metric flattening functions
	addpath('flatten_metric/')

	% Get the string for this text box
	target = get(hObject,'String');
	
    %extract name of source
    source_name = get(hObject,'Tag');
    
    %extract the column number as text
    source_number_text = source_name(end);
	
	% Check if this is a file, adding .mat on if needed
	if exist(target,'file') 
		
		filecheck = 1;
		
	elseif exist([target '.mat'],'file');
		
		filecheck = 1;
		original_target = target;
		target = [target '.mat'];
		
		
	else
		
		filecheck = 0;
		
	end
		
	% Different behavior based on filecheck
	if filecheck
		
		% if the file is a .mat file
		if strcmp('.mat',target(end-3:end))

			% Check if the file contains a convert structure
			target_contents = whos('-file',target);

			if strcmp('convert',{target_contents.name})

				% Enable the stretch button
				set(handles.(['stretch' source_number_text]),'Enable','on')

				% Set path to black text on white background
				set(hObject,'BackgroundColor','w','ForegroundColor','k')

			else

				% Disable the stretch button
				set(handles.(['stretch' source_number_text]),'Enable','off')


				% Set path to red text on white background to alert of lack of
				% convert structure
				set(hObject,'BackgroundColor','w','ForegroundColor','r')

			end
			
		else
			
			% Disable the stretch button
			set(handles.(['stretch' source_number_text]),'Enable','off')

			% The target is a non-.mat file
			if ~isempty(target)

				set(hObject,'BackgroundColor',[235 14 30]/255,'ForegroundColor','w')

			end	
			
		end
		
	else
		
		% Disable the stretch button
		set(handles.(['stretch' source_number_text]),'Enable','off')
		
		% If the target is nonempty, the path is not a file, so turn the
		% background red with black text
		if ~isempty(target)
		
			set(hObject,'BackgroundColor',[235 14 30]/255,'ForegroundColor','k')
	
		end
		
	end
	




end