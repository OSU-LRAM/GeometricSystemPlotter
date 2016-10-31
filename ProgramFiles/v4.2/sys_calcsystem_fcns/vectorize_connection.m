%Vectorize the connection definition functions
function s = vectorize_connection(s)

    %list of all components of the local connection that may be present
    connection_components = {'A_num','A_den','B_ref_point',...
        'B_rot_add','B_rot_mult'};
    
    %iterate over list of possible components
    for i = 1:length(connection_components)
        
        %check if component is present for this system
        if isfield(s,connection_components{i})
            
            %check if component is an inline function
            if isa(s.(connection_components{i}),'inline')
                
                %if it's inline, then vectorize it
                s.(connection_components{i}) = ...
                    vectorize(s.(connection_components{i}));
                
            %otherwise, check if the component is a function handle
            elseif isa(s.(connection_components{i}),'function_handle')
                
                %could put a test here to make sure the function gives
                %valid output
                
            else
                
                error(['The connection component ' connection_components{i}...
                    ' is neither an inline function nor a function handle'])
                
            end
            
        end
        
    end


end