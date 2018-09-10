% --- Enables or disables a block of checkboxes, based on an over-checkbox.
function enable_disable_block_Callback(source, eventdata, handles)
% source    handle of source checkbox
% handles   structure serving as name-to-handle lookup table
    
    %extract name of source
    source_name = get(source,'Tag');
    
    %text to insert into the button string
    switch source_name(1:2)
		case 'xy'
			inner_text = {'traj','net','BVI','cBVI'};
            switch source_name(3:5)
                case 'opt'
                    tailword = 'optcheckbox';
                otherwise
                    tailword = 'checkbox';
            end
           
		otherwise
			inner_text = {'X','Y','T','Xopt','Yopt','Topt'};
            tailword = 'checkbox';
	end
	
	
    %find beginning of string 'checkbox' in the source name
    insert_point = regexp(source_name,tailword);
    
    %loop over the text to be inserted, inserting it, then flipping the
    %enable bit on the resulting checkbox name
    for i = 1:length(inner_text)
        
        %insert the relevant text
        target_box = [source_name(1:insert_point-1) inner_text{i} source_name(insert_point:end)];
        
		%only try to enable or disable if the target box exists
		if isfield(handles,target_box)
		
			%enable or disable
			if get(source,'Value') && strcmp('on',get(source,'Enable'))

				set(handles.(target_box),'Enable','on');

			else

				set(handles.(target_box),'Enable','off');

			end
        
		end

	end
    
end