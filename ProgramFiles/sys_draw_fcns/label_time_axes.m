function label_time_axes(plot_info, component_list, ylabel_list)

	%Formatting options
	format_list = {'FontName','Helvetica','FontSize',24};

	%ylabel all plots
	for i = 1:length(plot_info.axes)

		%get which height function to use
		component_number = strmatch(plot_info.components{i}, component_list,'exact');

		%label the axis
		ylabel(plot_info.axes(i),ylabel_list{component_number},format_list{:},'Interpreter','latex');

	end



	%label the x values
	if plot_info.own_figure
		xlabel(plot_info.axes(end),'Time (s)',format_list{:});
	end

end