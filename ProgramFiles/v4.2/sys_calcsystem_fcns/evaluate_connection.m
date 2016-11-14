%Evaluate the connection over the fine grid for calculations and the coarse
%grid for vector display
function s = evaluate_connection(s)

    
    %list of zoom levels at which to evaluate connection vector fields and the
    %grid which should be used for them
    vector_field_list = {'display','vector';
                         'eval','eval'};
                     
    %list of all components of the local connection and metric that may be present
    connection_components = {'A_num','A_den','B_ref_point',...
        'B_rot_add','B_rot_mult','metric','metric_den'};
    
    
    %loop over list of zoom levels, creating the vector fields
    for i = 1:size(vector_field_list,1);
        
        % Shape values in grid
        a = s.grid.(vector_field_list{i,2});
        
        
        
        %iterate over list of possible components
        for j = 1:length(connection_components)

            %check if component is present for this system
            if isfield(s,connection_components{j})


                %generate one large array with all the vector information in it for
                %each of the numerator and denominator, and treat the set as two
                %cells in an array
                %note that the connection fields are the _negative_ of the
                %connection, this is handled in the merge_connection
                %function

                % Either call the connection function directly for the shape
                % inputs, or run a parallel loop to evaluate it at each of the
                % points
                switch s.function_type.(connection_components{j})

                    case 'block'
                        
                        % Evaluate the function at all points
                        A_full = s.(connection_components{j})(a{:});
                        
                        % Split the function output into pieces each the
                        % size of the input
                        A_cell_block = mat2tiles(A_full, size(a{1}));
                    
                    case 'woven'

                        % Evaluate the function at all points
                        A_full = s.(connection_components{j})(a{:});
                        
                        % Split the function output into pieces each the
                        % size of the output
                        a_point = cellfun(@(ai) ai(1),a,'UniformOutput',false);
                        A_cell = mat2tiles(A_full, size(s.(connection_components{j})(a_point{:})));
                        

                        % Swap the inner and outer structure of the cell so
                        % that outer layer is position in the tensor and
                        % inner layer is position in the grid.
                        A_cell_block = celltensorconvert(A_cell);

                    case 'single point'


                        A_cell = cell(size(a{1})); % Build a cell array to hold the function at each point
                        current_component = s.(connection_components{j});
                        parfor par_idx = 1:numel(a{1});   % Loop over all elements of the grid

                            % Extract the idx'th element of each grid, and put them
                            % into a cell array
                            a_point = cellfun(@(ai) ai(par_idx),a,'UniformOutput',false);

                            % Evaluate the function at the idx'th point
                            A_cell{par_idx} = current_component(a_point{:});
                        end
                        
                        % Swap the inner and outer structure of the cell so
                        % that outer layer is position in the tensor and
                        % inner layer is position in the grid.
                        A_cell_block = celltensorconvert(A_cell);


                    otherwise

                        error('Function type was specified as something other than block, woven, or single point')

                end
                       
                
                
                % Save the cell block evaluation of the function into the
                % appropriate field of the data structure
                
                if ~strncmp(connection_components{j},'metric',6)
                    
                    s.vecfield.(vector_field_list{i,1}).content.(connection_components{j})...
                        = A_cell_block;
                    
                    %mark what grid was used to create the field
                    s.vecfield.(vector_field_list{i,1}).type = (vector_field_list{i,2});
               else
                    
                    s.metricfield.(vector_field_list{i,1}).content.(connection_components{j})...
                        = A_cell_block;
                    
                    %mark what grid was used to create the field
                    s.metricfield.(vector_field_list{i,1}).type = (vector_field_list{i,2});
                    
                end
                
            end

        end
            
            
               
                
        
    end

end