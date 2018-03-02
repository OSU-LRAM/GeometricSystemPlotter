function plot_info = CCF_draw(s,p,plot_info,sys,shch,resolution)
%Draw the constraint curvature function
    
    %Get the configuration file, and extract the Colorpath
	configfile = 'sysplotter_config';
    configfile = fullfile(fileparts(mfilename('fullpath')),'..',configfile);
	load(configfile,'Colorset');

    %constraint curvature function list
    CCF_list = {'X','Y','T','Xopt','Yopt','Topt'};

    %Ensure that there are figure axes to plot into, and create new windows
    %for those axes if necessary
    plot_info = ensure_figure_axes(plot_info);  
	
	% Get the number of shape dimensions
	n_dim = numel(s.grid.eval);
	
	% get the number of position dimensions
	n_g = numel(s.DA);
    
	%if there is a singularity, deal with it for display
	if s.singularity
		
		CCF_addendum = '_scaled';
		
		singularity_location = logical(sum(cat(3,s.vecfield.eval.singularities{:}),3));
		%singularity_location = imdilate(singularity_location,[0 1 0;1 1 1;0 1 0]);
		
		%make the ztext reflect the arctangent scaling
		ztext = 'arctan(H)';
		
	else
		
		CCF_addendum = '';
		
		ztext = 'H';
		
	end

    
    %%%%%
    %constraint curvature function drawing
    
    %Extract the constraint curvature function
	
	%choose which constraint curvature function to plot
	switch plot_info.CCFtype
		
		case 'DA'
			
			H = cat(1,s.(['DA' CCF_addendum])...
				,s.(['DA_optimized' CCF_addendum]));
			
			if n_dim == 2
				
				h1name = 'DA';
				h2name = 'DA_optimized';
				
				h1 = cell(n_g,1);
				h2 = cell(n_g,1);
				
				if ~strcmp(shch,'null')
				
					for i = 1:numel(p.phi_locus)

						for j = 1:numel(p.phi_locus_full{i}.dA)

							h1{j}{i} = p.phi_locus_full{i}.(h1name){j};
							h2{j}{i} = p.phi_locus_full{i}.(h2name){j};

						end


					end

					zdata = cat(1,h1,h2);
				end
				
			end
			
		case 'dA'
			
			H = cat(1,s.(['dA' CCF_addendum])...
				,s.(['dA_optimized' CCF_addendum]));
			
			if n_dim == 2
				
				h1name = 'dA';
				h2name = 'dA_optimized';
				
				h1 = cell(n_g,1);
				h2 = cell(n_g,1);
				
				if ~strcmp(shch,'null')
					for i = 1:numel(p.phi_locus)

						for j = 1:numel(p.phi_locus_full{i}.dA)

							h1{j}{i} = p.phi_locus_full{i}.(h1name){j};
							h2{j}{i} = p.phi_locus_full{i}.(h2name){j};

						end

					end

					zdata = cat(1,h1,h2);
				end
				
			end			
			
			
		case 'LA'
						
			H = cat(1,cellfun(@(x,y) x-y,s.(['DA'  CCF_addendum])...
				,s.(['dA' CCF_addendum]),'UniformOutput',false)...
				,cellfun(@(x,y) x-y,s.(['DA_optimized' CCF_addendum])...
				,s.(['dA_optimized' CCF_addendum]),'UniformOutput',false));
			
			if n_dim == 2
				
				h1aname = 'DA';
				h1bname = 'dA';
				h2aname = 'DA_optimized';
				h2bname = 'dA_optimized';
				
				h1 = cell(n_g,1);
				h2 = cell(n_g,1);
				
				if ~strcmp(shch,'null')
					for i = 1:numel(p.phi_locus)

						for j = 1:numel(p.phi_locus_full{i}.dA)

							h1{j}{i} = p.phi_locus_full{i}.(h1aname){j}-p.phi_locus_full{i}.(h1bname){j};
							h2{j}{i} = p.phi_locus_full{i}.(h2aname){j}-p.phi_locus_full{i}.(h2bname){j};

						end

					end

					zdata = cat(1,h1,h2);
				end
				
			end						
		otherwise
			
			error('Unknown CCF function type')
			
	end
    
    %Extract the plotting grid
    grid = s.grid.eval;
    
	%%
	% Convert the function to the plotting grid specified in the gui
	[H,grid] = plotting_interp(H,grid,resolution,'scalar');
	
	%%%%%
	% If the shape coordinates should be transformed, make the conversion
	% (multiply the constraint curvature function by the inverse of the jacobian's
	% determinant)
	
	if plot_info.stretch == 1
			
		% Get the value by which to scale the constraint curvature function
		ascale = arrayfun(@(x,y) 1/det(s.convert.jacobian(x,y)),grid{:});

		% Apply the jacobian to the vectors
		H = cellfun(@(x) x.*ascale,H,'UniformOutput',false);

		% Convert the grid points to their new locations
		[grid{:}] = s.convert.old_to_new_points(grid{:});
		
	end
	
	
	


    %Make the plots
    for i = 1:length(plot_info.axes)
        
        %call up the relevant axis
        ax =plot_info.axes(i);
        
        %get which constraint curvature function to use
        function_number = strcmp(plot_info.components{i}, CCF_list);
        
		switch plot_info.style
			
			case 'surface'
		
% 				if s.singularity
% 					
% 					H{function_number}(singularity_location) = NaN;
% 					
% 				end
				
				%Plot the constraint curvature function
				meshhandle = surf(ax,grid{:},H{function_number});
				
				
				%If there's a shape change involved, plot it
				if ~strcmp(shch,'null')

					overlay_shape_change_3d_surf(ax,p,zdata{function_number,:},plot_info.stretch,s.convert,true);

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
				
				%Plot the constraint curvature function
				[junk, meshhandle] = contour(ax,grid{:},H{function_number},7,'linewidth',2);
				
				%If there's a shape change involved, plot it
				if ~strcmp(shch,'null')

					overlay_shape_change_2d(ax,p,plot_info.stretch,s.convert);

				end
				
				
				% Make edges if coordinates have changed
				if plot_info.stretch
					
					edgeres = 30;
					
					oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
						s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
					oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
						linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];

					[x_edge,y_edge] = s.convert.old_to_new_points(oldx_edge,oldy_edge);
					
					l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Color','k','LineWidth',1);
					
				end
				
				%Put an outline box around the plot
				box(ax,'on')

				%equal axes sized to match grid or new dimensions if
				%stretched
				
				if plot_info.stretch
					axis(ax,'equal');
					axis(ax,[min(grid{1}(:)) max(grid{1}(:)) min(grid{2}(:)) max(grid{2}(:))]);
				else
					axis(ax,'equal','tight');
				end
				
				%set the color map
				coloration = Colorset.colormap_contour;
                
                
            case 'pcolor'
                				
				%Plot the constraint curvature function
% 				[junk, meshhandle] = contour(ax,grid{:},H{function_number},7,'linewidth',2);

                

                [x_new, y_new] = s.convert.old_to_new_points(grid{1},grid{2});
                
                HF_isomap = interpn(s.convert.EI.A,s.convert.EI.B,s.convert.EI.C,x_new, y_new);
                
                isomap.x_new = x_new;
                isomap.y_new = y_new;
                isomap.HF_isomap = HF_isomap;
				
				%Plot the constraint curvature function
% 				[junk, meshhandle] = contour(ax,grid{:},H{function_number},7,'linewidth',2);
                meshhandle = pcolor(ax,grid{1},grid{2},H{function_number});
                
                meshhandle.ZData = HF_isomap;
                meshhandle.XData = x_new;
                meshhandle.YData = y_new;
                
                [grid{:}] = s.convert.old_to_new_points(grid{:});
				
				%If there's a shape change involved, plot it
				if ~strcmp(shch,'null')

					overlay_shape_change_2d(ax,p,plot_info.stretch,s.convert,isomap,s);

				end
				
				
				% Make edges if coordinates have changed
				if plot_info.stretch
					
					edgeres = 30;
					
					oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
						s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
					oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
						linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];

					[x_edge,y_edge] = s.convert.old_to_new_points(oldx_edge,oldy_edge);
                    
                    HF_isomap_edge = interpn(s.convert.EI.A,s.convert.EI.B,s.convert.EI.C,x_edge,y_edge);
					
					l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'ZData',HF_isomap_edge,'Color','k','LineWidth',1);
                    
%                     plot3(x_edge,y_edge,HF_isomap_edge,'k','Parent',ax,'LineWidth',1)
					
				end
				
				%Put an outline box around the plot
				box(ax,'on')

				%equal axes sized to match grid or new dimensions if
				%stretched
				
				if plot_info.stretch
					axis(ax,'equal');
					axis(ax,[min(grid{1}(:)) max(grid{1}(:)) min(grid{2}(:)) max(grid{2}(:))]);
				else
					axis(ax,'equal','tight');
				end
				
				%set the color map
				coloration = Colorset.colormap_contour;

				
			otherwise
				
				error('Unknown plot style for the constraint curvature function')
				
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
		label_shapespace_axes(ax,[],plot_info.stretch);

		%Set the tic marks
		set_tics_shapespace(ax,s,s.convert);



% 		%make hidden lines visible
% 		hidden2('off',ax)

		%%%%
		%Make clicking on the thumbnail open it larger in a new window

		if ~plot_info.own_figure

			%build a plot_info structure for just the current plot
			plot_info_specific.axes = 'new';
			plot_info_specific.components = plot_info.components(i);
			plot_info_specific.category = 'CCF';
			plot_info_specific.style = plot_info.style;
			plot_info_specific.CCFtype = plot_info.CCFtype;
			plot_info_specific.stretch = plot_info.stretch;

            %set the button down callback on the plot to be sys_draw with
			%the argument list for the current plot, and set the button
			%down callback for the mesh to the same
			set(plot_info.axes(i),'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});
			set(meshhandle,'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});

		else

			set(get(ax,'Parent'),'Name',[CCF_list{function_number} ' Constraint Curvature Function'])

			%Mark this figure as a constraint curvature function
			udata = get(plot_info.figure(i),'UserData');
			
			switch plot_info.style
				
				case 'surface'
					
					udata.plottype = 'CCF-surface';
					
				case 'contour'
					
					udata.plottype = 'CCF-contour';
					
			end
			
			set(plot_info.figure(i),'UserData',udata);
			
		end


        
    end
    
end