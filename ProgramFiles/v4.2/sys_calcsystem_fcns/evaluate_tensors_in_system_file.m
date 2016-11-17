function s = evaluate_tensors_in_system_file(s,component_list,zoom_list,destination)

    %loop over list of zoom levels, creating the vector fields
    for i = 1:size(zoom_list,1);
        
        % Shape values in grid
        a = s.grid.(zoom_list{i,2});
        
        
        
        %iterate over list of possible components
        for j = 1:length(component_list)

            %check if component is present for this system
            if isfield(s,component_list{j})


                %Evaluate the function over the grid, as a cell array
                %where the outer structure is that of the tensor, and the
                %inner structure is that of the grid
                
                s.(destination).(zoom_list{i,1}).content.(component_list{j}) =...
                    evaluate_tensor_over_grid(s.(component_list{j}),a);
                
                
                % Mark what zoom level was used to create this field
                s.(destination).(zoom_list{i,1}).type = zoom_list{i,2};
                
            end

        end
            
            
               
                
        
    end
    
end
