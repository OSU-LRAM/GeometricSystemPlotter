function plotpublishpushbutton_enable_menu(handles)
%Enable or disable the plot pushbutton, based on the values of the system
%and shape change menus

    %Get the system and shape change info
    system_index = get(handles.systemmenu,'Value');
    system_names = get(handles.systemmenu,'UserData');
    
    current_system = system_names{system_index};
    
    shch_index = get(handles.shapechangemenu,'Value');
    shch_names = get(handles.shapechangemenu,'UserData');
    
    current_shch = shch_names{shch_index};
    
    %get a list of all plot and publish pushbuttons
    all_elements = fieldnames(handles);
    plot_ind = strmatch('plotpushbutton',all_elements);
    publish_ind = strmatch('publishpushbutton',all_elements);
    plot_or_publish_ind = [plot_ind;publish_ind];
    plot_or_publish = all_elements(plot_or_publish_ind);

    %If either system or shch is a default, disable the plot and publish
    %buttons
    if strcmp(current_system,'default') || strcmp(current_shch,'default')
        
        button_target_state = 'off';
    
    %otherwise, enable the buttons
    else
        
        button_target_state = 'on';
        
    end
    
    for i = 1:length(plot_or_publish)

        set(handles.(plot_or_publish{i}),'Enable',button_target_state);

    end
    
end