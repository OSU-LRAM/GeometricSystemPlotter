% visualizes predicted third-order bounds
% designed to augment existing gait plots
function [len, cBVI_fun, to_fun] = add_bound_plots(s, shape, P, phi, do_opt, fig_start)
    % sanitize
    if ~exist('phi', 'var')
        phi = 0;
    end
    if ~exist('do_opt', 'var')
        do_opt = true;
    end
    if ~exist('fig_start', 'var')
        fig_start = 1;
    end
    if ~ishandle(fig_start)
        error('Open trajectory graph (or gait graph, in addition) before adding third-order bounds');
    end
    
    % compute bound
    [len, cBVI_fun, to_fun] = bound_third_order(s, shape,...
                                                @(a,b) A_est_center(a,b,do_opt),...
                                                @(a,b) cBVI_est_taylor(a,b,do_opt),...
                                                P, phi);    
    % infer what to plot based on the number of open figures
    if ishandle(fig_start+1)
        draw_gait_bound(fig_start, shape, len);
        draw_traj_bound(fig_start+1, P, cBVI_fun, to_fun, [0 len]);
    else
        draw_traj_bound(fig_start, P, cBVI_fun, to_fun, [0 len]);
    end
end

% draws a circle, centered at <shape>, of diameter l
function draw_gait_bound(fig, shape, l)
    % construct vertices
    t = linspace(0, 2*pi);
    points = shape(:) + l/2*[cos(t); sin(t)];
    circle_poly = polyshape(points(1,:), points(2,:));
    % plotting
    figure(fig);
    hold on
    poly = plot(circle_poly);
    poly.EdgeColor = [1 1 1]*0.4;
    poly.FaceColor = [1 1 1]*0.4;
    poly.FaceAlpha = 0.5;
    poly.LineStyle = '--';
    hold off
end

% draws third-order/cBVI prediction over <range> (use s.grid_range)
% TODO: rethink how this is done; almost certainly not right
function draw_traj_bound(fig, P, cBVI_fun, to_fun, range)
    % construct vertices
    l_vec = linspace(range(1), range(2));
    to_points = cell2mat(arrayfun(to_fun, l_vec, 'UniformOutput', false));
    cBVI_points = -cBVI_fun(l_vec);
    plot_points = [cBVI_points + to_points, flip(cBVI_points - to_points,2)];
    shape = polyshape(plot_points(1,:), plot_points(2,:));
    % plotting
    figure(fig);
    hold on
    poly = plot(shape, 'DisplayName', ['3^{rd} Order Bound, p=', num2str(P)]);
    poly.EdgeColor = [1 1 1]*0.4;
    poly.FaceColor =[1 1 1]*0.4;
    poly.FaceAlpha = 0.5;
    poly.LineStyle = '--';
    hold off
end

