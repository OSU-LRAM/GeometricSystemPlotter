%Evaluate the connection over the fine grid for calculations and the coarse
%grid for vector display
function s = evaluate_metric_derivatives(s)

    % First, extract the metric and the grid
    M = s.metricfield.metric_eval.content.metric;
    grid = s.grid.metric_eval;
    baseline = grid_to_baseline(grid);
    
    % Now construct an nx1 cell array to hold the derivatives of M with
    % respect to the shape variables
    dM = repmat({M},numel(grid),1);
    
    %%%%%%
    % Now fill the cell array with the gradient of M (using ndgradient
    % because it matches with ndgrid)
    
    % Cell array to temporarily hold the output of gradient
    grad_components = cell(numel(grid),1);
    
    % Loop over every element of the metric
    for idx_metric = 1:numel(M)
        
        % Get the gradient of that component, using ndmetric to correspond
        % with the ndgrid used by our code
        [grad_components{:}] = ndgradient(M{idx_metric},baseline{:});
        
        % Loop over all derivative components, inserting the gradient
        % component into the appropriate location
        for idx_component = 1:numel(dM)
            dM{idx_component}{idx_metric} = grad_components{idx_component};
        end
        
    end
    
    %%%%%%%%%%%%%%%%
    % Second derivative
    
    % Now construct an nxn cell array to hold the derivatives of M with
    % respect to the shape variables
    ddM = repmat({M},numel(grid));
    
    % Loop over every element of the metric
    for idx_metric = 1:numel(M)
        
       % Loop over every element of the gradient
       for idx_component1 = 1:numel(dM)
           
           % Get the gradient of that component of the gradient
           [grad_components{:}] = ndgradient(dM{idx_component1}{idx_metric},baseline{:});
           
           % Loop over all components of this derivative, inserting the
           % gradient into the appropriate location
           for idx_component2 = 1:numel(dM)
               ddM{idx_component1,idx_component2}{idx_metric} = grad_components{idx_component2};
           end
       end
       
    end

    %%%%%%%%%%%%%%%%
    % Finally, save the output into the structure
    s.coriolisfield.coriolis_eval.content.dM = dM;
    s.coriolisfield.coriolis_eval.content.ddM = ddM;
    
    
%    dmdalpha = shape_partial_mass(M_full,J_full,local_inertias,jointangles,A_eval,A_grid)
    
%     %list of all components of the local connection and metric that may be present
%     component_list = {'metric','metric_den'};
%     
%     
%     % Evaluate all components in the list
%     
%     s = evaluate_tensors_in_system_file(s,component_list,{'metric_eval','metric_eval'},'metricfield');
%     
%     % resample metric at a resolution appropriate for displaying as an
%     % ellipse field
%     s = resample_tensors_in_system_file(s,component_list,'metric_eval',{'metric_display','metric_display'},'metricfield');
    

end