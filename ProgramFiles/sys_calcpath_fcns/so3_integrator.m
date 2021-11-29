% so3_integrator
% computes displacement over time for a kinematic system on SO(3)
% inputs:
    % t_vec: vector of timesteps
    % alpha_fun: vector-valued function of time, producing shape variables
    % d_alpha_fun: time derivative of alpha
    % grid: ndgrid for A
    % A: local connection tensor
% outputs:
    % g_circ: exponential coordinates of displacements
    % g: cell of xyz transform matrices with respect to timesteps t

function [g_circ, g] = so3_integrator(t_vec, alpha_fun, d_alpha_fun, grid, A)
    % cayley fn for orthogonalization
    cay  = @(A) (eye(3)+A)\(eye(3)-A);
    % preallocate cell of matrices
    g = cell(1,length(t_vec));
    g{1} = eye(3);
    % for each timestep (gap in t):
    for i = 2:length(t_vec)
        % get shape alpha, ss velocity d_alpha
        t = t_vec(i);
        alpha_vec = alpha_fun(t);
        alpha_cell = num2cell(alpha_vec);
        d_alpha_vec = d_alpha_fun(t);
        % get LC matrix A at alpha (interpolate)
        A_mat = zeros(size(A));
        for entry = 1:numel(A)
            A_mat(entry) = interpn(grid{:},A{entry},alpha_cell{:});
        end
        % use d_alpha to produce xyz velocity
        g_circ_mat = vec_to_mat_SO3(A_mat * d_alpha_vec);
        % exponentiate
        g_mat_t = expm((t_vec(i)-t_vec(i-1))*g_circ_mat);
        % cumprod with previous cell
        g{i} = g{i-1}*g_mat_t;
        % re-orthogonalize
        %g{i} = g{i}/det(g{i}); %bad method
        g{i} = cay(cay(g{i}));
    end
    % take matrix log of ea. cell for return
    g_circ = zeros(3,length(t_vec));
    for i = 1:length(t_vec)
        g_circ(:,i) = mat_to_vec_SO3(logm(g{i}));
    end
end