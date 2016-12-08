function [x_min,x_max,y_min,y_max] = ...
	set_axis_limits(target_axis,xdata,ydata,xbuffer,ybuffer)
	%set the axis limits for the plot defined by xdata and ydata

	%set the x range
	x_min = min(xdata(:)) - (xbuffer * range(xdata(:)));
	x_max = max(xdata(:)) + (xbuffer * range(xdata(:)));

	%set the y range
	y_min = min(ydata(:)) - (ybuffer * range(ydata(:)));
	y_max = max(ydata(:)) + (ybuffer * range(ydata(:)));

	%if either in fact has no range, set a default range of +/- 1
	if abs(x_max-x_min) < 10*eps

		x_center = x_min;

		x_min = x_center-1;
		x_max = x_center+1;

	end
	if abs(y_max-y_min) < 10*eps

		y_center = y_min;

		y_min = y_center-1;
		y_max = y_center+1;

	end

	%apply the axis limits
	axis(target_axis,[x_min,x_max,y_min,y_max]);


end