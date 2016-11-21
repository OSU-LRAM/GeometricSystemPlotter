function s = calc_stretch_functions(s)
% Apply cartographic normalization to the shape space

    s.convert = fast_flatten_metric(s.grid.metric_display,s.metricfield.metric_display.content.metric);
    
end