function plot_info = ShapeSpace_draw(s,p,plot_info,sys,shch,resolution)

    if strcmp(plot_info.components{1},'opt')
        s.geometry.baseframe = sys;
    else
        % do nothing
    end
    
    plot_info = ensure_figure_axes(plot_info)
    

    illustrate_shapespace(s,plot_info.axes)

end
