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