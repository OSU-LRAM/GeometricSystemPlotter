% estimates total lie bracket body vel vector using grid value at center
function cBVI_est = cBVI_est_center(s, shape, do_opt)
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
    % interpolate grid to evaluate DA at desired shape
    DA_vec = zeros(3,1); %assume SE(2)
    for dim = 1:3
        DA_vec(dim) = interp2(s.grid.eval{2}, s.grid.eval{1},...
                              DA{dim},...
                              shape(1), shape(2));
    end
    cBVI_est = @(l) pi/4 * l^2 * DA_vec;
end