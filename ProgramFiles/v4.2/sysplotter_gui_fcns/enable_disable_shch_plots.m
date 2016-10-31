%Disable shape change dependent functions if No shape change is
%selected, and enable them if there is a shape change selected (or there is
%no selection)
function enable_disable_shch_plots(source,eventdata,handles)


    
    %get shape change menu state
    shch_index = get(handles.shapechangemenu,'Value');
    shch_names = get(handles.shapechangemenu,'UserData');
    
    current_shch = shch_names{shch_index};
    
    %get names of handles
    handlenames = fieldnames(handles);
    
    %If "No shape change" is selected, then disable the shape change
    %dependent plots
    if strcmp(current_shch,'null')
        
        shch_dependent_enable = 'off';
        
    else
        
        shch_dependent_enable = 'on';
        
    end
        
    %classes of plot dependent on shape changes
    plot_types = {'bvi','disp','xy','xyopt'};
    
    %get all top level checkbox handles for the shape change dependent
    %plots
    top_handles = [];
    for i = 1:length(plot_types)
    
        top_handles = [top_handles; strmatch([plot_types{i} 'checkbox'],handlenames)];
        
    end
    
    %loop over the handles, setting the enable field on the corresponding
    %checkboxes, and then calling the block enable for their children
    for i = 1:length(top_handles)
        
        %set the enable field
        set(handles.(handlenames{top_handles(i)}),'Enable',shch_dependent_enable);
        
        %call the block handler
        enable_disable_block_Callback(handles.(handlenames{top_handles(i)}), eventdata, handles);
        
    end
    
end