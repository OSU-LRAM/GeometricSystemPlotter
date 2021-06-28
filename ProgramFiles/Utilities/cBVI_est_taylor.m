% estimates total lie bracket body vel vector using taylor approximation
function cBVI_est = cBVI_est_taylor(s, shape, do_opt)
    % sanitize
    if ~exist('do_opt', 'var')
        do_opt = true;
    end
    % get total lie bracket
    if do_opt
        DA = s.DA_optimized;
    else
        DA = s.DA;
    end
    % for each dimension of DA, compute gradient 
    [f_x, f_y, f_xx, f_xy, f_yy, f_yx] = deal(cell(size(DA)));
    for dim = 1:3
        baseline = grid_to_baseline(s.grid.eval);
        [f_x{dim}, f_y{dim}] = ndgradient(DA{dim}, baseline{:});
        [f_xx{dim}, f_xy{dim}] = ndgradient(f_x{dim}, baseline{:});
        baseline_rev = flip(baseline);
        [f_yx{dim}, f_yy{dim}] = ndgradient(f_y{dim}, baseline_rev{:});
    end
    % get taylor ingredients (yes, I'm aware that this is awful)
    [f_i, f_x_i, f_y_i, f_xx_i, f_xy_i, f_yy_i, f_yx_i] = deal(zeros(3,1));
    for dim = 1:3
        f_i(dim) = grad_interp(s, DA{dim}, shape);
        f_x_i(dim) = grad_interp(s, f_x{dim}, shape);
        f_y_i(dim) = grad_interp(s, f_y{dim}, shape);
        f_xx_i(dim) = grad_interp(s, f_xx{dim}, shape);
        f_xy_i(dim) = grad_interp(s, f_xy{dim}, shape);
        f_yx_i(dim) = grad_interp(s, f_yx{dim}, shape);
        f_yy_i(dim) = grad_interp(s, f_yy{dim}, shape);
    end
    % generate taylor poly
    T = @(dx,dy) f_i +... %fn value
                 f_x_i*dx + f_y_i*dy +... %first order
                 1/2*(f_xx_i*dx^2 + (f_xy_i+f_yx_i)*dx*dy + f_yy_i*dy^2);
    % produce expression w.r.t. diameter
    syms dr dphi
    T_p = T(dr*cos(dphi), dr*sin(dphi));
    T_r = int(T_p, dphi, [0 2*pi]);
    cBVI_r = matlabFunction(int(T_r*dr, dr)); %radius-dependent cBVI
    cBVI_est = @(len) cBVI_r(len/2); %cast to diameter
end

function scalar = grad_interp(s, f_dim, x)
    scalar = interp2(s.grid.eval{2}, s.grid.eval{1}, f_dim, x(1), x(2));
end