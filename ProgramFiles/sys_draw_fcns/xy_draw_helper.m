function plot_info = xy_draw_helper(s,p,plot_info,sys,shch,optselect);
% Draw the displacement locus indicated by optselect
	
	%Ensure that there are figure axes to plot into, and create new windows
	%for those axes if necessary
	plot_info = ensure_figure_axes(plot_info);
	
    %Get the configuration file, and extract the Colorpath
	configfile = './sysplotter_config';
	load(configfile,'Colorset');
	
	% define a color list
	colorlist = {Colorset.spot, Colorset.secondary, 'k'};
	
	% Get the axes to plot into
	ax = plot_info.axes;
	
	% Prime arrays for net x and y displacement
	net_x = zeros(numel(p.G_locus_full),1);
	net_y = zeros(numel(p.G_locus_full),1);
	
	% Extract list of components to include
	components = plot_info.components{1};
	
	% prime the line element handles
	[traj_h,BVI_h,cBVI_h,net_h] = deal([]);
	
	% prime the accumulators
	[traj_x,BVI_x,cBVI_x,net_x,traj_y,BVI_y,cBVI_y,net_y] = deal([]);
	
	% prime the legend strings
	[traj_l,BVI_l,cBVI_l,net_l] = deal({});
	
	%%%%
	% Loop over all gaits
	for i = 1:numel(p.G_locus_full)
		
		% Set the line color
		if  numel(p.G_locus_full) == 1
			traj_color = Colorset.spot;
			BVI_color = Colorset.spot;
			cBVI_color = Colorset.spot;
		else
			listnum = mod(i-1,numel(colorlist))+1;
			[traj_color,BVI_color,cBVI_color] = deal(colorlist{listnum});
		end
		
		% draw the trajectory
		if any(strncmp('traj',components,4))
			traj_h(i) = line('Parent',ax,'XData',p.G_locus_full{i}.(['G' optselect])(:,1)...
				,'YData',p.G_locus_full{i}.(['G' optselect])(:,2),'LineWidth',3,'Color',traj_color);
		
			% Store the data to use in axis-limit calculation
			traj_x = [traj_x; p.G_locus_full{i}.(['G' optselect])(:,1)];
			traj_y = [traj_y; p.G_locus_full{i}.(['G' optselect])(:,2)];
			
			% provide a legend string
			if  numel(p.G_locus_full) == 1
				lstring = 'Trajectory';
			else
				lstring = ['Trajectory ' num2str(i)];
			end
			
			traj_l = [traj_l;lstring];
			
		end
		

		% draw the bvi point
		if any(strncmp('BVI',components,3))
			BVI_h(i) = line('Parent',ax,'XData',p.G_locus_full{i}.(['bvi' optselect])(end,1)...
				,'YData',p.G_locus_full{i}.(['bvi' optselect])(end,2)...
				,'Marker','o','MarkerSize',16,'LineWidth',3,'Color',BVI_color,'Linestyle','none');

			BVI_x = [BVI_x; p.G_locus_full{i}.(['bvi' optselect])(end,1)];
			BVI_y = [BVI_y; p.G_locus_full{i}.(['bvi' optselect])(end,2)];
			
			
			% provide a legend string
			if  numel(p.G_locus_full) == 1
				lstring = 'BVI';
			else
				lstring = ['BVI ' num2str(i)];
			end
			
			BVI_l = [BVI_l;lstring];
		end
		
		% draw the cbvi point
		if any(strncmp('cBVI',components,4)) 
			if isfield(p,['cBVI' optselect]) && ~isempty(p.(['cBVI' optselect]){i})
				cBVIx = p.(['cBVI' optselect]){i}(1);
				cBVIy = p.(['cBVI' optselect]){i}(2);
			else
				cBVIx = [];
				cBVIy = [];
			end
			cBVI_h(i) = line('Parent',ax,'XData',cBVIx...
				,'YData',cBVIy...
				,'Marker','x','MarkerSize',16,'LineWidth',3,'Color',cBVI_color,'Linestyle','none');
		
			cBVI_x = [cBVI_x; cBVIx];
			cBVI_y = [cBVI_y; cBVIy];
			
			% provide a legend string
			if  (numel(p.G_locus_full) == 1) 
				lstring = 'cBVI';
			else
				lstring = ['cBVI ' num2str(i)];
			end
			
			if (~(isfield(p,['cBVI' optselect]) && ~isempty(p.(['cBVI' optselect])))) || (numel(s.grid.vector) > 2)
				lstring = [lstring ' (not calculated)'];
			end
			
			cBVI_l = [cBVI_l;lstring];

		end		
		
		net_x(i) = p.G_locus_full{i}.(['G' optselect])(end,1);
		net_y(i) = p.G_locus_full{i}.(['G' optselect])(end,2);
		
		
	end
	
	
	% draw the net displacement locus for all the gaits
	if any(strncmp('net',components,3))
		net_color = 'k';
		net_h = line('Parent',ax,'XData',net_x...
			,'YData',net_y,'ZData',ones(size(net_x)),'LineWidth',3 ...
			,'Marker','o','MarkerFaceColor',net_color,'MarkerSize',12,'Color',net_color);
		
		net_l = {'Net displacement'};
	end
	
	%%%%%%%%
	% Set the axis limits
	
	if any(~isempty([traj_h(:); BVI_h(:); cBVI_h(:); net_h(:)]))
		% Collect all the displacements
		x_collect = cat(1,traj_x(:),BVI_x(:),cBVI_x(:),net_x(:),0);
		y_collect = cat(1,traj_y(:),BVI_y(:),cBVI_y(:),net_y(:),0);

		
		
		% Set the axes with a nice buffering
		[x_min,x_max,y_min,y_max] = ...
			set_axis_limits(ax,x_collect,y_collect,.07,.07);
    end
    set(ax,'ZLim',10*[-1 1]); % prevent layering from going outside zlim
	axis(ax,'equal');
	new_lim = [get(ax,'xlim') get(ax,'ylim')];
    
	
	
	%%%%%%%
	% Draw the axis lines
	xaxis_h(1) = line('Parent',ax','YData',[0 0],'XData',[0 new_lim(2)],'ZData',-2*[1 1]...
		,'LineWidth',2,'LineStyle','--','Color','k');

	xaxis_h(2) = line('Parent',ax','YData',[0 0],'XData',[0 new_lim(1)],'ZData',-2*[1 1]...
		,'LineWidth',2,'LineStyle','--','Color','k');
	
	yaxis_h(1) = line('Parent',ax','XData',[0 0],'YData',[0 new_lim(4)],'ZData',-2*[1 1]...
		,'LineWidth',2,'LineStyle','--','Color','k');

	yaxis_h(2) = line('Parent',ax','XData',[0 0],'YData',[0 new_lim(3)],'ZData',-2*[1 1]...
		,'LineWidth',2,'LineStyle','--','Color','k');
	
	
	%%%%%%%
	% Turn on the outline box
	box(ax,'on')
    
	
	%%%%%
	% Insert a legend
	L = legend([traj_h,BVI_h,cBVI_h,net_h],traj_l{:},BVI_l{:},cBVI_l{:},net_l{:});

	%%%%
	%Make clicking on the thumbnail open it larger in a new window

	if ~plot_info.own_figure

		%build a plot_info structure for just the current plot
		plot_info_specific.axes = 'new';
		plot_info_specific.components = plot_info.components;
		if isempty(optselect)
			plot_info_specific.category = 'xy';
		else
			plot_info_specific.category = 'xyopt';
		end
		plot_info_specific.style = plot_info.style;

		%set the button down callback on the plot to be sys_draw with
		%the argument list for the current plot, and set the button
		%down callback for the mesh to the same
		set(plot_info.axes,'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});

	else
		
		%Mark this figure as an xy locus plot
		udata = get(plot_info.figure,'UserData');
		
		if isempty(optselect)
			set(get(ax(1),'Parent'),'Name','XY original locus')
			udata.plottype = 'xy';
		else
			set(get(ax(1),'Parent'),'Name','XY optimal locus')
			udata.plottype = 'xyopt';
		end
		
		set(plot_info.figure,'UserData',udata,'Renderer','painters');

	end
end