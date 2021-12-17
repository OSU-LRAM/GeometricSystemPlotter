function s = test_connection_and_metric(s)

%list of all components of the local connection and metric that may be present
    connection_components = {'A_num','A_den','B_ref_point',...
        'B_rot_add','B_rot_mult','metric','metric_den'};
    
    %iterate over list of possible components
    for i = 1:length(connection_components)
        
        %check if component is present for this system
        if isfield(s,connection_components{i})
            
            
            %%%%%%
            % Verify that the connection numerator is a function
            
            % check if the component is a function handle
            if isa(s.(connection_components{i}),'function_handle')
                
                %could put a test here to make sure the function gives
                %valid output

            % check if component is an inline function
            elseif isa(s.(connection_components{i}),'inline')
                
                %if it's an inline function, then vectorize it
                s.(connection_components{i}) = ...
                    vectorize(s.(connection_components{i}));
                
                
            else
                
                error(['The connection component ' connection_components{i}...
                    ' is neither an inline function nor a function handle'])
                
           
            end
            
            
            %%%%%%%%
            % Test the functions to find out what kind of output they use
            s = test_function_type(s,connection_components{i});
            
            
        end
        
    end


end