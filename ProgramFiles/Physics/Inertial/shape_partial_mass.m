function dmdalpha = shape_partial_mass(geometry,physics,jointangles,A_eval,A_grid)
[~, ~, ~, J_full, ~, M_full, local_inertias] = Inertial_connection_discrete(geometry,physics,jointangles);
dJdq = mobile_jacobian_derivative(J_full);
dMdq = partial_mass_matrix(J_full,dJdq,local_inertias,'mobile');
dmdalpha = pull_back_partial_mass(M_full,dMdq,A_eval,A_grid,jointangles);
end


function dmdalpha = pull_back_partial_mass(M,dMdq,A,A_grid,jointangles)
    angle_cells = num2cell(jointangles);
    num_joints = size(A,2);
    dMdalpha = dMdq(4:end); % Get the partials for the joints, skipping the first three terms
                            % (Since these are the g-circ terms)
    % Get the local connection at the given joint angles
    A_interp = cellfun(@(C) interpn(A_grid{:},C,jointangles(1),jointangles(2),'spline'),A);
    % Get dA/dalpha at the given jointangles
    dA_interp = cell(num_joints,1);
    % TODO: Need some way to make this number-of-joints-agnostic. Note
    % that results are given as [dA/dalpha2 dA/dalpha1] because the
    % first term is difference in horizontal direction and second is
    % difference in vertical direction; the grid has alpha1 vary by
    % column and alpha2 vary by row.
    [dA2, dA1] = cellfun(@(C) gradient(C,A_grid{2}(1,:),A_grid{1}(:,1)),A,'UniformOutput',false);
    dA_interp{1} = cellfun(@(C) interpn(A_grid{:},C,angle_cells{:},'spline'),dA1);
    dA_interp{2} = cellfun(@(C) interpn(A_grid{:},C,angle_cells{:},'spline'),dA2);
    
    dmdalpha = cell(num_joints,1);
    for i = 1:numel(dmdalpha)
        
        dmdalpha{i} = [-dA_interp{i}', zeros(num_joints)]*M*[-A_interp; eye(num_joints)] ...
            + [-A_interp', eye(num_joints)]*dMdalpha{i}*[-A_interp; eye(num_joints)] ...
            + [-A_interp', eye(num_joints)]*M*[-dA_interp{i}; zeros(num_joints)];
    end
end