function plot_info = hfun_draw(s,p,plot_info,sys,shch,convert,resolution)
%Draw the height function
    
    %Get the configuration file, and extract the Colorpath
	configfile = './sysplotter_config';
	load(configfile,'Colorset');

    %height function list
    hfun_list = {'X','Y','T','Xopt','Yopt','Topt'};

    %Ensure that there are figure axes to plot into, and create new windows
    %for those axes if necessary
    plot_info = ensure_figure_axes(plot_info);  
	
	% Get the number of shape dimensions
	n_dim = numel(s.grid.eval);
	
	% get the number of position dimensions
	n_g = numel(s.height);
    
	%if there is a singularity, deal with it for display
	if s.singularity
		
		hfun_addendum = '_scaled';
		
		singularity_location = logical(sum(cat(3,s.vecfield.eval.singularities{:}),3));
		%singularity_location = imdilate(singularity_location,[0 1 0;1 1 1;0 1 0]);
		
		%make the ztext reflect the arctangent scaling
		ztext = 'arctan(H)';
		
	else
		
		hfun_addendum = '';
		
		ztext = 'H';
		
	end

    
    %%%%%
    %Height function drawing
    
    %Extract the height function
	
	%choose which height function to plot
	switch plot_info.hfuntype
		
		case 'cH'
			
			H = cat(1,s.(['height_corrected' hfun_addendum])...
				,s.(['height_optimized_corrected' hfun_addendum]));
			
			if n_dim == 2
				
				h1name = 'height_corrected';
				h2name = 'height_optimized_corrected';
				
				h1 = cell(n_g,1);
				h2 = cell(n_g,1);
				
				if ~strcmp(shch,'null')
				
					for i = 1:numel(p.phi_locus)

						for j = 1:numel(p.phi_locus_full{i}.height)

							h1{j}{i} = p.phi_locus_full{i}.(h1name){j};
							h2{j}{i} = p.phi_locus_full{i}.(h2name){j};

						end


					end

					zdata = cat(1,h1,h2);
				end
				
			end
			
		case 'oH'
			
			H = cat(1,s.(['height' hfun_addendum])...
				,s.(['height_optimized' hfun_addendum]));
			
			if n_dim == 2
				
				h1name = 'height';
				h2name = 'height_optimized';
				
				h1 = cell(n_g,1);
				h2 = cell(n_g,1);
				
				if ~strcmp(shch,'null')
					for i = 1:numel(p.phi_locus)

						for j = 1:numel(p.phi_locus_full{i}.height)

							h1{j}{i} = p.phi_locus_full{i}.(h1name){j};
							h2{j}{i} = p.phi_locus_full{i}.(h2name){j};

						end

					end

					zdata = cat(1,h1,h2);
				end
				
			end			
			
			
		case 'dH'
						
			H = cat(1,cellfun(@(x,y) x-y,s.(['height_corrected'  hfun_addendum])...
				,s.(['height' hfun_addendum]),'UniformOutput',false)...
				,cellfun(@(x,y) x-y,s.(['height_optimized_corrected' hfun_addendum])...
				,s.(['height_optimized' hfun_addendum]),'UniformOutput',false));
			
			if n_dim == 2
				
				h1aname = 'height_corrected';
				h1bname = 'height';
				h2aname = 'height_optimized_corrected';
				h2bname = 'height_optimized';
				
				h1 = cell(n_g,1);
				h2 = cell(n_g,1);
				
				if ~strcmp(shch,'null')
					for i = 1:numel(p.phi_locus)

						for j = 1:numel(p.phi_locus_full{i}.height)

							h1{j}{i} = p.phi_locus_full{i}.(h1aname){j}-p.phi_locus_full{i}.(h1bname){j};
							h2{j}{i} = p.phi_locus_full{i}.(h2aname){j}-p.phi_locus_full{i}.(h2bname){j};

						end

					end

					zdata = cat(1,h1,h2);
				end
				
			end						
		otherwise
			
			error('Unknown height function type')
			
	end
    
    %Extract the plotting grid
    grid = s.grid.eval;
    
	%%
	% Convert the function to the plotting grid specified in the gui
	[H,grid] = plotting_interp(H,grid,resolution,'scalar');
	
	%%%%%
	% If the shape coordinates should be transformed, make the conversion
	% (multiply the height function by the inverse of the jacobian's
	% determinant)
	
	if ~isempty(convert)
			
		% Get the value by which to scale the height function
		ascale = arrayfun(@(x,y) 1/det(convert.jacobian(x,y)),grid{:});

		% Apply the jacobian to the vectors
		H = cellfun(@(x) x.*ascale,H,'UniformOutput',false);

		% Convert the grid points to their new locations
		[grid{:}] = convert.old_to_new_points(grid{:});
		
	end
	
	
	


    %Make the plots
    for i = 1:length(plot_info.axes)
        
        %call up the relevant axis
        ax =plot_info.axes(i);
        
        %get which height function to use
        function_number = strmatch(plot_info.components{i}, hfun_list,'exact');
        
		switch plot_info.style
			
			case 'surface'
		
% 				if s.singularity
% 					
% 					H{function_number}(singularity_location) = NaN;
% 					
% 				end
				
				%Plot the height function
				meshhandle = surf(ax,grid{:},H{function_number});
				
				
				%If there's a shape change involved, plot it
				if ~strcmp(shch,'null')

					overlay_shape_change_3d_surf(ax,p,zdata{function_number,:},convert);

				end

				%square axes
				axis(ax,'tight')
				

				%Put an outline box around the plot
				nicebox(ax,'on')
				
				%load the color map
                coloration = Colorset.colormap;
				if s.singularity
					coloration = coloration.^5;
				end

				
			case 'contour'
				
				%Plot the height function
				[junk, meshhandle] = contour(ax,grid{:},H{function_number},7,'linewidth',2);
				
				%If there's a shape change involved, plot it
				if ~strcmp(shch,'null')

					overlay_shape_change_2d(ax,p,convert);

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
				
				%Put an outline box around the plot
				box(ax,'on')

				%equal axes sized to match grid or new dimensions if
				%stretched
				
				if isempty(convert)
					axis(ax,'equal');
					axis(ax,[min(grid{1}(:)) max(grid{1}(:)) min(grid{2}(:)) max(grid{2}(:))]);
				else
					axis('equal','tight');
				end
				
				%set the color map
				coloration = Colorset.colormap_contour;

				
			otherwise
				
				error('Unknown plot style for the height function')
				
		end

		%Iterate up the tree to find the figure that contains the current
		%axis
		parH = get(ax,'Parent');
		while 1

			if strmatch('figure',get(parH,'Type'))

				break

			else

				parH = get(parH,'Parent');

			end

		end

		set(parH,'Colormap',coloration);

		%center the color map around zero
		Clim = get(ax,'Clim'); %get the current color limits
		C_outer = max(abs(Clim)); %get the maximum distance from zero
		set(ax,'Clim',[-C_outer C_outer]); %set an inclusive range around zero



		%Label the axes
		label_shapespace_axes(ax,[],~isempty(convert));

		%Set the tic marks
		set_tics_shapespace(ax,s,convert);



% 		%make hidden lines visible
% 		hidden2('off',ax)

		%%%%
		%Make clicking on the thumbnail open it larger in a new window

		if ~plot_info.own_figure

			%build a plot_info structure for just the current plot
			plot_info_specific.axes = 'new';
			plot_info_specific.components = plot_info.components(i);
			plot_info_specific.category = 'hfun';
			plot_info_specific.style = plot_info.style;
			plot_info_specific.hfuntype = plot_info.hfuntype;
			plot_info_specific.stretch = plot_info.stretch;
			plot_info_specific.stretchpath = plot_info.stretchpath;

			%set the button down callback on the plot to be sys_draw with
			%the argument list for the current plot, and set the button
			%down callback for the mesh to the same
			set(plot_info.axes(i),'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch});
			set(meshhandle,'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch});

		else

			set(get(ax,'Parent'),'Name',[hfun_list{function_number} ' Height Function'])

			%Mark this figure as a height function
			udata = get(plot_info.figure(i),'UserData');
			
			switch plot_info.style
				
				case 'surface'
					
					udata.plottype = 'hfun-surface';
					
				case 'contour'
					
					udata.plottype = 'hfun-contour';
					
			end
			
			set(plot_info.figure(i),'UserData',udata);
			
		end


        
    end
    
end