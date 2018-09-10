function plot_info = xy_draw(s,p,plot_info,sys,shch,resolution)
% Standin for drawing the unoptimized xy plot

    if regexp(plot_info.components{1}{1},'opt$')
        optselect = '_opt';
    else
        optselect = '';
    end
    
	plot_info = xy_draw_helper(s,p,plot_info,sys,shch,optselect);

end