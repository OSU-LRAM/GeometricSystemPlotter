function dmdalpha = shape_partial_mass(M_full,J_full,local_inertias,jointangles,A_eval,A_grid)
% Calculates the partial of the pulled-back mass matrix with respect to the
% shape variables alpha_i. This is done by calculating the partial of the
% full mass matrix with respect to the full set of configuration variables
% q = [g_circ, r]', where g_circ represents the world configuration and r
% represents the local shape variables.
%
% Inputs:
%   geometry: Structure containing information about geometry of chain. See
%       Inertial_connection_discrete for information about required fields.
%
%   physics: Structure containing information about the system's physics.
%       See Inertial_connection_discrete for information about required 
%       fields.
%
%   jointangles: Array containing values of shape positions at which the
%       pulled-back partial mass matrix should be evaluated
%
%   A_eval: Cell array containing the local connection evaluated at
%       discrete points determined by the values in A_grid. As of 1/21/20, 
%       this is given by s.vecfield.eval.content.A_num
%
%   A_grid: Points at which the local connection in A_eval was evaluated.
%       As of 1/21/20, this is given by s.grid.eval
%
% Output for a system with m links and k joints:
%
%   dmdalpha: Pulled-back partial of the mass matrix, in terms of the shape
%       variables only. (k by 1) cell array where the ith entry corresponds
%       to the partial of the pulled-back mass matrix with respect to the
%       ith shape variable.

% Get the full jacobian and mass matrices as well as the links' local
% inertia tensors

% Obtain the derivative of the Jacobian with respect to the configuration
dJdq = mobile_jacobian_derivative(J_full);
% Get the partial of the mass matrix with respect to the configuration
% variables
dMdq = partial_mass_matrix(J_full,dJdq,local_inertias);
% Pull back the partial of the mass matrix to be in terms of the shape
% variables only
dmdalpha = pull_back_partial_mass(M_full,dMdq,A_eval,A_grid,jointangles);
end


function dmdalpha = pull_back_partial_mass(M,dMdq,A,A_grid,jointangles)
% Pulls back the full partial mass matrix dMdq to be in terms of the shape
% variables only.
%
% Inputs for a system with m links and q of size (n by 1):
%   M: Full mass matrix of size (n by n)
%
%   dMdq: Partial of full mass matrix with respect to the configuration
%       variables q; (n by 1) cell array where the ith cell contains an (n by n)
%       matrix corresponding to dMdq_i, or the partial of the mass matrix with
%       respect to the ith configuration variable in q
%
%   A_eval: Cell array containing the local connection evaluated at
%       discrete points determined by the values in A_grid. As of 1/21/20, 
%       this is given by s.vecfield.eval.content.A_num
%
%   A_grid: Points at which the local connection in A_eval was evaluated.
%       As of 1/21/20, this is given by s.grid.eval
%
%   jointangles: Array containing values of shape positions at which the
%       pulled-back partial mass matrix should be evaluated
%
% Output for a system with m links and k joints:
%
%   dmdalpha: Pulled-back partial of the mass matrix, in terms of the shape
%       variables only. (k by 1) cell array where the ith entry corresponds
%       to the partial of the pulled-back mass matrix with respect to the
%       ith shape variable.

    angle_cells = num2cell(jointangles);
    num_joints = size(A,2);
    dMdalpha = dMdq(4:end); % Get the partials for the joints, skipping the first three terms
                            % (Since these are the g-circ terms)
    % Get the local connection at the given joint angles
    ja = num2cell(jointangles);
    A_interp = cellfun(@(C) interpn(A_grid{:},C,ja{:},'spline'),A);
    % Get dA/dalpha at the given jointangles
    dA_interp = cell(num_joints,1);
    % TODO: Need some way to make this number-of-joints-agnostic. Note
    % that results are given as [dA/dalpha2 dA/dalpha1] because the
    % first term is difference in horizontal direction and second is
    % difference in vertical direction; the grid has alpha1 vary by
    % column and alpha2 vary by row.
    
    A_grid_mesh = A_grid([2,1,3:numel(A_grid)]);
    
    grid_baseline = cell(size(A_grid_mesh));
    callout_template = num2cell(ones(size(grid_baseline)));
    for idx = 1:numel(grid_baseline)
        callout = callout_template;
        callout{idx} = ':';
        grid_baseline{idx} = squeeze(A_grid{idx}(callout{:}));
    end

    dA = repmat({cell(size(A))},1,numel(jointangles));
    [dA{[2,1,3:numel(A_grid)]}] = cellfun(@(C) gradient(C,grid_baseline{[2,1,3:numel(A_grid)]}),A,'UniformOutput',false);
    
    for idx = 1:numel(dA)
        dA_interp{idx} = cellfun(@(C) interpn(A_grid{:},C,angle_cells{:},'spline'),dA{idx});
    end
%     dA_interp{2} = cellfun(@(C) interpn(A_grid{:},C,angle_cells{:},'spline'),dA2);

%     [dA2, dA1] = cellfun(@(C) gradient(C,A_grid{2}(1,:),A_grid{1}(:,1)),A,'UniformOutput',false);
%     dA_interp{1} = cellfun(@(C) interpn(A_grid{:},C,angle_cells{:},'spline'),dA1);
%     dA_interp{2} = cellfun(@(C) interpn(A_grid{:},C,angle_cells{:},'spline'),dA2);
    
    dmdalpha = cell(num_joints,1);
    % The pulled-back mass matrix m_alpha is given by the matrix product:
    %       m_alpha = [-A' I]*M*[A; I]
    % So the partial of the pulled-back mass matrix is found simply by
    % performing the chain rule on this quantity, as is done below.
    for i = 1:numel(dmdalpha)
        
        dmdalpha{i} = [-dA_interp{i}', zeros(num_joints)]*M*[-A_interp; eye(num_joints)] ...
            + [-A_interp', eye(num_joints)]*dMdalpha{i}*[-A_interp; eye(num_joints)] ...
            + [-A_interp', eye(num_joints)]*M*[-dA_interp{i}; zeros(num_joints)];
    end
end