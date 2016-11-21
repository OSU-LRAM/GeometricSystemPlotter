function plot_info = xyopt_draw(s,p,plot_info,sys,shch,resolution)
% Standin for drawing the unoptimized xy plot

	plot_info = xy_draw_helper(s,p,plot_info,sys,shch,'_opt');

end