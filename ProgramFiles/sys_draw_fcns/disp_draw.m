function plot_info = disp_draw(s,p,plot_info,sys,shch,convert,resolution) %#ok<INUSD>
%draw the bvi plots

    %Get the configuration file, and extract the Colorpath
	configfile = './sysplotter_config';
	load(configfile,'Colorset');
	
	% define a color/linestyle list
	linelist = {Colorset.spot, '-';
		Colorset.secondary, '--';
		'k' ,':';
		Colorset.spot, '--';
		Colorset.secondary, ':';
		'k' ,'-'};
	
% 	switch plot_info.linestyle
% 		
% 		case 'mono'
% 			
% 			linelist = {Colorset.spot, '-'};
% 			
% 		case 'cycle'
% 	
%             linelist = {Colorset.spot, '-';
%                 Colorset.secondary, '--';
%                 'k' ,':';
%                 Colorset.spot, '--';
%                 Colorset.secondary, ':';
%                 'k' ,'-'};
% 		
% 		otherwise
% 			
% 			error('Unexpected value for plot_info.linestyle')			
% 			
% 	end

	% Number of position space dimensions
	n_dim = size(p.G_locus{1}{1}.G,2);

	%BVI function list (matching a list in the gui code)
	disp_list = {'X','Y','T','Xopt','Yopt','Topt'};

	%Ensure that there are figure axes to plot into, and create new window
	%for those axes if necessary
	plot_info = ensure_figure_subaxes(plot_info);

	%%%
	%Draw the displacement

	%make the plots
	for i = 1:length(plot_info.axes)

		%call up the relevant axis
		ax = plot_info.axes(i);

		%make the plot boxed in from all sides
		box(ax,'on')

		%get which bvi to use
		disp_number = mod(strmatch(plot_info.components{i}, disp_list,'exact')-1,n_dim)+1;
		
		% Decide which component to plot
		disp_field = ['G_' plot_info.components{i}(2:end)];
		if strcmp(disp_field(end),'_')
			disp_field = disp_field(1:(end-1));
		end
		
		% loop over all paths
		for j = 1:numel(p.G_locus_full)
		
			% select a linestyle from the list
			listnum = mod(j-1,size(linelist,1))+1;

			% Plot the BVI lines
			line('Parent',ax,'XData',p.time_full{j}...
				,'YData',p.G_locus_full{j}.(disp_field)(:,disp_number)...
				,'LineWidth',2.5,'Color',linelist{listnum,1},'LineStyle',linelist{listnum,2});
			
		end
		
		%%%%
		%set the limits for pretty display
		% extract all the y data
		ydata = cellfun(@(x) x.(disp_field)(:,disp_number),p.G_locus_full,'UniformOutput',false);
		[x_min,x_max,y_min,y_max] = ...
			set_axis_limits(plot_info.axes(i),cat(1,p.time_full{:})...
			,cat(1,ydata{:}),0.0,0.07);

		%put in the divisions -- the ends of the time segments for the
		%first path
		divide_time = cellfun(@(x) max(x),p.time{1});
		vertical_dividers(ax,divide_time,y_min,y_max);

		%Set the tic marks
		set_tics_time_history(ax,s,(i == length(plot_info.axes)));
		
		
		%%%%
		%Make clicking on the thumbnail open it larger in a new window

		if ~plot_info.own_figure

			%build a plot_info structure for just the current plot
			plot_info_specific.axes = 'new';
			plot_info_specific.components = plot_info.components;
			plot_info_specific.category = 'disp';
			plot_info_specific.linestyle = plot_info.linestyle;

			%if this window is not an optimized-system window, strip out
			%any optimized plot calls, otherwise, split out unoptimized
			%calls
			if sum(strcmp(plot_info.components{i},disp_list(1:3)))

				split_ind = find( (strcmp(disp_list{1},plot_info.components)) | ...
					(strcmp(disp_list{2},plot_info.components)) | ...
					(strcmp(disp_list{3},plot_info.components)));

			else

				split_ind = find( (strcmp(disp_list{4},plot_info.components)) | ...
					(strcmp(disp_list{5},plot_info.components)) | ...
					(strcmp(disp_list{6},plot_info.components)));

			end

			if ~isempty(split_ind)
				plot_info_specific.components = plot_info_specific.components(split_ind);
			end

			%set the button down callback on the plot to be sys_draw with
			%the argument list for the current plot, and set the button
			%down callback for the mesh to the same
			set(ax,'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});

		else

			set(get(ax(1),'Parent'),'Name','Displacements')

			%Mark this figure as a Displacement
			udata = get(plot_info.figure,'UserData');
			udata.plottype = 'disp';
			set(plot_info.figure,'UserData',udata);


		end

	end

	%label the axes
ylabel_list = {'$x$','$y$','$\theta$','$x$','$y$','$\theta$'};	label_time_axes(plot_info, disp_list,ylabel_list);


end