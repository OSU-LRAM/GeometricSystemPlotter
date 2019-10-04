%Evaluate inertial properties over the fine grid for calculations and the coarse
%grid for vector display
function s = evaluate_inertial_properties(s)

    %list of all components of the local connection and metric that may be present
    component_list = {{'M_alpha'},{'dM_alphadalpha'}};
    eval_list = {{'mass_eval','mass_eval'},{'coriolis_eval','coriolis_eval'}};
    display_list = {{'mass_display','mass_display'},{'coriolis_display','coriolis_display'}};
    field_list = {'massfield','coriolisfield'};
    
    % Evaluate all components in the list
    for i = 1:length(component_list)
        s = evaluate_tensors_in_system_file(s,component_list{i},eval_list{i},field_list{i});

%         % resample metric at a resolution appropriate for displaying as an
%         % ellipse field
%         s = resample_tensors_in_system_file(s,component_list{i},eval_list{i}{1},display_list{i},field_list{i});
    end
    % Use a numerical approach for the second partial derivative of the
    % mass matrix
    s = evaluate_mass_second_derivative_numerical(s);
end

function s = evaluate_mass_second_derivative_numerical(s)
    % This function assumes that the first partial derivative
    % dM_alpha/dalpha has already been calculated over the grid specified
    % by s.grid.coriolis_eval and stored in s.coriolisfield. Note that,
    % because this is the second partial derivative, each partial
    % derivative will need to have its own partial derivative taken with
    % respect to each of the shape variables, which augments the space of
    % the second partial derivative to be one greater than that of the
    % partial derivative.
    
    % Square cell with dimension equal to the size of the dM_alphadalpha
    % vector, since dimension is augmented
    num_joints = length(s.coriolisfield.coriolis_eval.content.dM_alphadalpha);
    grid = s.grid.coriolis_eval;
    ddM_alphadalpha = cell(num_joints);
    for first_partial = 1:num_joints
        dM_first_partial = s.coriolisfield.coriolis_eval.content.dM_alphadalpha{first_partial};
        [dA2, dA1] = cellfun(@(C) gradient(C,grid{2}(1,:),grid{1}(:,1)),dM_first_partial,'UniformOutput',false);
        ddM_alphadalpha{first_partial,1} = dA1;
        ddM_alphadalpha{first_partial,2} = dA2;
    end
    s.coriolisfield.coriolis_eval.content.ddM_alphadalpha = ddM_alphadalpha;
end