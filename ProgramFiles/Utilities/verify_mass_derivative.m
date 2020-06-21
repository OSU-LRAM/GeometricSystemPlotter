function valid = verify_mass_derivative(s,err_tol)
    % Set a default error tolerance if one wasn't provided
    if ~exist('err_tol','var') || isempty(err_tol)
        err_tol = 1e-6;
    end
    % Pick a random shape in the range of valid shapes
    b = s.grid_range(2);
    a = s.grid_range(1);
    shape = (b-a).*rand(1,2) + a;
    shapelist = num2cell(shape);
    
    % Evaluate M at the random shape
    M = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    
    % Evaluate M at a small delta_alpha in both alpha_1 and alpha_2
    % directions to verify if dM_alphadalpha is accurate
    dalpha = 0.001;
    shapelist_a1 = shapelist;
    shapelist_a1{1} = shapelist_a1{1} + dalpha;
    M_a1 = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist_a1{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    shapelist_a2 = shapelist;
    shapelist_a2{2} = shapelist_a2{2} + dalpha;
    M_a2 = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist_a2{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    
    % Get what dM_alphadalpha says from analytical precalculation
    dM_alphadalpha = cell(size(shapelist));
    for i = 1:length(shapelist)
        dM_alphadalpha{i} = cellfun(@(C) interpn(s.grid.coriolis_eval{:},C,...
            shapelist{:},'spline'),s.coriolisfield.coriolis_eval.content.dM_alphadalpha{i});
    end
    
    % Check the difference between dM_alphadalpha and the numerical deltas
    err_a1 = (M_a1-M)./dalpha - dM_alphadalpha{1};
    err_a2 = (M_a2-M)./dalpha - dM_alphadalpha{2};
    
    if any(any(abs(err_a1) > err_tol)) || any(any(abs(err_a2) > err_tol))
        valid = false;
    else
        valid = true;
    end
end