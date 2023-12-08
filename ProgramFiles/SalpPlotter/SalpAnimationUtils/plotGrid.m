function plotGrid(ax,base_anim_geom)
% Plot grid marks to distinguish salp motion

    ticks = base_anim_geom.ticks;

    for idx = 1:numel(ticks.points{1})
        x_ticks = ticks.points{1};
        plot(ax,[x_ticks(idx),x_ticks(idx)], ...
            [ticks.plot_min{2},ticks.plot_max{2}], ...
            'Color', ticks.color);
    end

    for idx = 1:numel(ticks.points{2})
        y_ticks = ticks.points{2};
        plot(ax,[ticks.plot_min{1},ticks.plot_max{1}], ...
            [y_ticks(idx),y_ticks(idx)], ...
            'Color', ticks.color);
    end
end