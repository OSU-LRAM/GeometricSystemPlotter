% optimize_so3_sys.m
% sysplotter-relevant wrapper for general optimize_so3 function

function s = optimize_so3_sys(s)
    % pull out grid spacing as cell of vectors
    grid = s.grid.finite_element;
    samples = cellfun(@(x) unique(x), grid, 'UniformOutput', false);
    % optimize
    A_orig = s.vecfield.eval.content.Avec;
    [~, A_opt, X, Y, Z, grad_X, grad_Y, grad_Z] = optimize_so3(samples, A_orig);
    % save results in sys struct
    s.vecfield.eval.content.Avec_optimized = A_opt;
    s.B_optimized.eval.Beta = {X, Y, Z};
    s.B_optimized.eval.gradBeta = [grad_X; grad_Y; grad_Z];
end