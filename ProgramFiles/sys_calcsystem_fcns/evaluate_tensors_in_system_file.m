function s = evaluate_tensors_in_system_file(s,component_list,zoom_list,destination)

    %loop over list of zoom levels, creating the vector fields
    for i = 1:size(zoom_list,1)
        
        % Shape values in grid
        a = s.grid.(zoom_list{i,2});
        
        
        
        %iterate over list of possible components
        for j = 1:length(component_list)

            %check if component is present for this system
            if isfield(s,component_list{j})

                % Check to see if the user has acknowledged the presence of
                % a singularity in the connection and or metric
                ignore_singularity_warning = isfield(s,'ignore_singularity_warning') && s.ignore_singularity_warning;
                

                %Evaluate the function over the grid, as a cell array
                %where the outer structure is that of the tensor, and the
                %inner structure is that of the grid
                
               
                   T =  evaluate_tensor_over_grid(s.(component_list{j}),a,ignore_singularity_warning);
                
                % If there are any NaN values and a singularity has not been called out
                % inpaint the NaNs and warn the user

                 nan_present = cellfun(@(Tc) any(isnan(Tc(:))),T);
                if ~isfield(s,'singularity') && ...
                    any(nan_present(:))
                    

                    T = cellfun(@(Tc) inpaint_nans(Tc,4),T,'UniformOutput',false);
                    
                    if ~ignore_singularity_warning
                        warning('NaN values were inpainted on a tensor, but a singularity was not specified in the file')
                    end
                end

                 s.(destination).(zoom_list{i,1}).content.(component_list{j}) = T;
                
                % Mark what zoom level was used to create this field
                s.(destination).(zoom_list{i,1}).type = zoom_list{i,2};
                                
            end

        end
        
    end
    
end