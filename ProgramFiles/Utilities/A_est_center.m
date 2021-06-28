% estimates local connection matrix using grid value at center
function A_mat = A_est_center(s, shape, do_opt)
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
    % interpolate grid to evaluate LC at desired shape
    A_mat = zeros(size(A));
    for dim = 1:3
        for shvar = 1:length(shape)
            contrib = interp2(s.grid.eval{2}, s.grid.eval{1},...
                              A{dim, shvar},...
                              shape(1), shape(2));
            A_mat(dim, shvar) = contrib;
        end
    end
end