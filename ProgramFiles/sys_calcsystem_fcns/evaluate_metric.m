%Evaluate the connection over the fine grid for calculations and the coarse
%grid for vector display
function s = evaluate_metric(s)

    
    %list of all components of the local connection and metric that may be present
    component_list = {'metric','metric_den'};
    
    
    % Evaluate all components in the list
    
    s = evaluate_tensors_in_system_file(s,component_list,{'metric_eval','metric_eval'},'metricfield');
    
    % resample metric at a resolution appropriate for displaying as an
    % ellipse field
    s = resample_tensors_in_system_file(s,component_list,'metric_eval',{'metric_display','metric_display'},'metricfield');
    

end