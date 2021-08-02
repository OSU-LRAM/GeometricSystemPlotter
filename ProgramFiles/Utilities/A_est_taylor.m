function A_est = A_est_taylor(s, shape, do_opt)
    % sanitize
    if ~exist('do_opt', 'var')
        do_opt = true;
    end
    % get local connection
    if do_opt
        A = s.vecfield.eval.content.Avec_optimized;
    else
        A = s.vecfield.eval.content.Avec;
    end
    
    % we're going to generate a taylor polynomial for every LC vecfield.
    % effectively we map a gait radius to a 3x2 LC matrix
    
    % for each vector field, compute gradient(s) w.r.t. iension(s)
    [f_x, f_y, f_xx, f_yy, f_xy, f_yx] = deal(cell(size(A)));
    for i = 1:numel(A)
        baseline = grid_to_baseline(s.grid.eval);
        [f_x{i}, f_y{i}] = ndgradient(A{i}, baseline{:});
        [f_xx{i}, f_xy{i}] = ndgradient(f_x{i}, baseline{:});
        baseline_rev = flip(baseline);
        [f_yx{i}, f_yy{i}] = ndgradient(f_y{i}, baseline_rev{:});
    end
    
    % interpolate gradients at desired shape (this is still bad, but fast)
    [f_i,f_x_i,f_y_i,f_xx_i,f_xy_i,f_yy_i,f_yx_i] = deal(zeros(size(A)));
    for i = 1:numel(A)
        f_i(i) = grad_interp(s, A{i}, shape);
        f_x_i(i) = grad_interp(s, f_x{i}, shape);
        f_y_i(i) = grad_interp(s, f_y{i}, shape);
        f_xx_i(i) = grad_interp(s, f_xx{i}, shape);
        f_xy_i(i) = grad_interp(s, f_xy{i}, shape);
        f_yx_i(i) = grad_interp(s, f_yx{i}, shape);
        f_yy_i(i) = grad_interp(s, f_yy{i}, shape);
    end
    
    % construct taylor polynomial
    T = @(dx,dy) f_i +... %fn value
                 f_x_i*dx + f_y_i*dy +... %first order
                 1/2*(f_xx_i*dx^2 + (f_xy_i+f_yx_i)*dx*dy + f_yy_i*dy^2);
    % produce expression w.r.t. diameter
    syms dr dphi
    T_p = T(dr*cos(dphi), dr*sin(dphi)); % polar polynomial
    phi_samp = linspace(0, 2*pi, 3); % sample radius
    A_samp_full = subs(T_p, dphi, phi_samp);
    A_samp = reshape(A_samp_full, size(A,1), size(A,2), length(phi_samp));
    A_mean = mean(A_samp, 3);    
    A_r = matlabFunction(A_mean);
    
    A_est = @(len) A_r(len/2);
end

function scalar = grad_interp(s, f_i, x)
    scalar = interp2(s.grid.eval{2}, s.grid.eval{1}, f_i, x(1), x(2));
end