function plot_info = vfield_draw(s,p,plot_info,sys,shch,convert,resolution)

    %Vector field list
    vfield_list = {'X','Y','T','Xopt','Yopt','Topt'};
	
	% Get the number of dimensions
	n_dim = numel(s.grid.eval);
    
    %Ensure that there are figure axes to plot into, and create new windows
    %for those axes if necessary
    plot_info = ensure_figure_axes(plot_info);    
    
    %%%%%%%
    % Get the vector field and interpolate into the specified grid
	
    % Extract the display vector field
    V = cat(1,s.vecfield.display.content.Avec,s.vecfield.display.content.Avec_optimized);
    
    % Extract the plotting grid
    grid = s.grid.vector;
	
	%%
	% Convert the vector field to the plotting grid specified in the gui
	[V,grid] = plotting_interp(V,grid,resolution,'vector');
	
	%%%%%
	% If the shape coordinates should be transformed, make the conversion
	% (multiply the vectors by the inverse jacobian)
	
	if ~isempty(convert)
	
		% Calculate the jacobians at the plotting points
		Jac = arrayfun(convert.jacobian,grid{:},'UniformOutput',false);
		
		% Use the jacobians to convert the vectors
		
		% Iterate over all connection vector fields present
		for i = 1:size(V,1)
			
			% Iterate over all vectors present
			for j = 1:numel(V{i,1})
				
				% Extract all components of the relevant vector
				tempVin = cellfun(@(x) x(j),V(i,:));
				
				% Multiply by the inverse Jacobian
				tempVout = Jac{j}\tempVin(:);
				
				% Replace vector components
				for k = 1:size(V,2)
					
					V{i,k}(j) = tempVout(k);
					
				end
				
			end
			
		end

		% Convert the grid points to their new locations
		[grid{:}] = convert.old_to_new_points(grid{:});
		
	end
    
    %%%
    %If there's a singularity, use arctan scaling on the magnitude of the
    %vector field
    if s.singularity
    
		V = arctan_scale_vector_fields(V);
        
	end
    
	
	
    
    for i = 1:length(plot_info.axes)
        
        %call up the relevant axis
        ax = plot_info.axes(i);
        
        %get which vector field to use
        field_number = strmatch(plot_info.components{i}, vfield_list,'exact');
        
		%plot the vector field arrows
        if n_dim == 2
			quiver(ax,grid{:},V{field_number,1},V{field_number,2},'k','LineWidth',2)
		else
			quiver3(ax,grid{1:3},V{field_number,:},'k','LineWidth',2)
		end
			
		% Make edges if coordinates have changed
		if ~isempty(convert)

			edgeres = 30;

			oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
				s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
			oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
				linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];

			[x_edge,y_edge] = convert.old_to_new_points(oldx_edge,oldy_edge);

			l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Color','k','LineWidth',1);

		end
				
		%set the display range
		if isempty(convert)
			axis(ax,s.grid_range);
		end
		
		if isempty(convert)
			axis(ax,'equal');
			axis(ax,[min(grid{1}(:)) max(grid{1}(:)) min(grid{2}(:)) max(grid{2}(:))]);
		else
			axis('equal','tight');
		end
        
        %Label the axes (two-dimensional)
        label_shapespace_axes(ax,[],~isempty(convert));
		
        %Set the tic marks
        set_tics_shapespace(ax,s,convert);
        
        %If there's a shape change involved, plot it
        if ~strcmp(shch,'null')
						
            overlay_shape_change_2d(ax,p,convert);
            
        end
        
        
        %%%%
        %Make clicking on the thumbnail open it larger in a new window
        
        if ~plot_info.own_figure

            %build a plot_info structure for just the current plot
            plot_info_specific.axes = 'new';
            plot_info_specific.components = plot_info.components(i);
            plot_info_specific.category = 'vfield';
			plot_info_specific.stretch = plot_info.stretch;
			plot_info_specific.stretchpath = plot_info.stretchpath;
			
            %set the button down callback on the plot to be sys_draw with
            %the argument list for the current plot
            set(plot_info.axes(i),'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch});
            
        else
            
            set(get(ax(1),'Parent'),'Name',[vfield_list{field_number} ' Vector Field'])
			
			%Mark this figure as a vector field
			udata = get(plot_info.figure(i),'UserData');
			udata.plottype = 'vfield';
			set(plot_info.figure(i),'UserData',udata);

        
        end
        
        
    end
    
end