function formatSalpPlot(ax,salp,base_anim_geom)
 
    %Set axis bounds
    minX = base_anim_geom.ticks.plot_min{1};
    maxX = base_anim_geom.ticks.plot_max{1};
    minY = base_anim_geom.ticks.plot_min{2};
    maxY = base_anim_geom.ticks.plot_max{2};

    % Fit axis so we're always zoomed on the swimmer
    axis(ax,[minX,maxX,minY,maxY]);
    if isfield(base_anim_geom,'overrideAxisLimits')
        if base_anim_geom.overrideAxisLimits
            axis(ax,[minX,maxX,-.2,.2]);
        end
    end

    title(ax,salp.name);
    set(ax,'XTick',[0]);
    set(ax,'YTick',[0]);
end