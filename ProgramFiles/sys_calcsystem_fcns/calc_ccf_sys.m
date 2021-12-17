% calc_ccf_sys.m
% sysplotter wrapper for constraint curvature function

function s = calc_ccf_sys(s)
    % get grid points in ea. dim
    grid = s.grid.eval;
    samples = cellfun(@(x) unique(x), grid, 'UniformOutput', false);
    % get matrix, vector repr. fns
    mat_fn = s.conf_space.mat_fn;
    vec_fn = s.conf_space.vec_fn;
    % compute CCF for orig. opt coords
    A_orig = s.vecfield.eval.content.Avec;
    A_opt = s.vecfield.eval.content.Avec_optimized;
    [s.DA, s.dA, s.lb] = calc_ccf(samples, A_orig, mat_fn, vec_fn);
    [s.DA_optimized, s.dA_optimized, s.lb_optimized] = calc_ccf(samples, A_opt, mat_fn, vec_fn);
end