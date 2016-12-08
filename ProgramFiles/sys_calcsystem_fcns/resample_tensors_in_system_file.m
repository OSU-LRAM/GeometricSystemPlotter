function s = resample_tensors_in_system_file(s,component_list,base_zoom,zoom_list,destination)

    %loop over list of zoom levels, creating the vector fields
    for i = 1:size(zoom_list,1);
        
        % Shape values in grid
        a = s.grid.(zoom_list{i,2});
        
        
        
        %iterate over list of possible components
        for j = 1:length(component_list)

            %check if component is present for this system
            if isfield(s,component_list{j})

                % Create resampled field
                s.(destination).(zoom_list{i,1}).content.(component_list{j})...
                    = cellfun(@(V) interpn(s.grid.(base_zoom){:},V,a{:})...
                    ,s.(destination).(base_zoom).content.(component_list{j}),'UniformOutput',false);                
                
                % Mark what zoom level was used to create this field
                s.(destination).(zoom_list{i,1}).type = zoom_list{i,2};
                
            end

        end
            
            
               
                
        
    end
    
end