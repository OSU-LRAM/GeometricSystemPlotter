function grid_baseline = grid_to_baseline(grid)

    grid_baseline = cell(size(grid));
    callout_template = num2cell(ones(size(grid_baseline)));
    for idx = 1:numel(grid_baseline)
        callout = callout_template;
        callout{idx} = ':';
        grid_baseline{idx} = squeeze(grid{idx}(callout{:}));
    end
    
end