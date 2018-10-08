%Evaluate the connection over the fine grid for calculations and the coarse
%grid for vector display
function s = evaluate_connection(s)

    
                     
    %list of all components of the local connection and metric that may be present
    component_list = {'A_num','A_den','B_ref_point',...
        'B_rot_add','B_rot_mult'};
    
    
    % Evaluate all components in the list at the eval density
    
    s = evaluate_tensors_in_system_file(s,component_list,{'eval','eval'},'vecfield');
    
    % resample the fields at a low density for vector field display
    s = resample_tensors_in_system_file(s,component_list,'eval',{'display','vector';'finite_element','finite_element'},'vecfield');
        
    
    


end