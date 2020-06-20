function plot_info = beta_draw(s,p,plot_info,sys,shch,resolution)
%Draw the height function

    %Get the configuration file, and extract the Colorpath
    configfile = './sysplotter_config';
    load(configfile,'Colorset');

	%height function list
	beta_list = {'X','Y','T'};
	beta_names = {'$\beta_{x}$','$\beta_{y}$','$\beta_{\theta}$'};

	%Ensure that there are figure axes to plot into, and create new windows
	%for those axes if necessary
	plot_info = ensure_figure_axes(plot_info);

	%%%%%%
	% Extract the beta values as zdata
	
	% Get the number of shape dimensions
	n_dim = numel(s.grid.eval);
	
	% get the number of position dimensions
	n_g = numel(s.dA);
    
    if n_dim ==2 
        % Set up the zdata for the plot
        if ~strcmp(shch,'null')
            if n_dim == 2

                zdata = cell(n_g,1);

                for i = 1:numel(p.phi_locus)

                    for j = 1:numel(p.phi_locus_full{i}.Beta)

                        zdata{j}{i} = p.phi_locus_full{i}.Beta{j};

                    end

                end

            end
        end	

        %%%%%
        %Height function drawing

        %Extract the height function
        B = s.B_optimized.eval.Beta;

        %Extract the plotting grid
        grid = s.grid.eval;

        %%
        % Convert the function to the plotting grid specified in the gui
        [B,grid] = plotting_interp(B,grid,resolution,'scalar');

%         %%%%%
%         % No need to change the value of the beta function for stretching,
%         % because it is a zero-form, not a 2-form
% 
%         if plot_info.stretch && (numel(s.grid.eval) == 2)
% 
%             % Convert the grid points to their new locations
%             [grid{:}] = s.convert.old_to_new_points(grid{:});
% 
%         end



        %Make the plots
        for i = 1:length(plot_info.axes)

            %call up the relevant axis
            ax =plot_info.axes(i);

            %get which height function to use
            function_number = strcmp(plot_info.components{i}, beta_list);

            switch plot_info.style

                case 'surface'

                    %Plot the height function
                    if n_dim==2
                        meshhandle = surf(ax,grid{:},B{function_number});
                    end

                    %If there's a shape change involved, plot it
                    if ~strcmp(shch,'null')

                        overlay_shape_change_3d_surf(ax,p,zdata{function_number,:},plot_info.stretch,s.convert,false);

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
                    [junk, meshhandle] = contour(ax,grid{:},B{function_number},7,'linewidth',2);

                    %If there's a shape change involved, plot it
                    if ~strcmp(shch,'null')

                        overlay_shape_change_2d(ax,p,plot_info.stretch,s.convert);

                    end


%                     % Make edges if coordinates have changed
%                     if plot_info.stretch && (numel(s.grid.eval) == 2)
% 
%                         edgeres = 30;
% 
%                         oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
%                             s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
%                         oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
%                             linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];
% 
%                         [x_edge,y_edge] = s.convert.old_to_new_points(oldx_edge,oldy_edge);
% 
%                         l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Color','k','LineWidth',1);
% 
%                     end

                    %Put an outline box around the plot
                    box(ax,'on')

                    %square axes
                    axis(ax,'equal','tight')

                    %set the color map
                    coloration = Colorset.colormap_contour;

                otherwise

                    error('Unknown plot style for the beta function')

            end

            %Iterate up the tree to find the figure that contains the current
            %axis
            parB = get(ax,'Parent');
            while 1

                if strmatch('figure',get(parB,'Type'))

                    break

                else

                    parB = get(parB,'Parent');

                end

            end
            set(parB,'Colormap',coloration);

            %center the color map around zero
            Clim = get(ax,'Clim'); %get the current color limits
            C_outer = max(abs(Clim)); %get the maximum distance from zero
            set(ax,'Clim',[-C_outer C_outer]); %set an inclusive range around zero


            %square axes
            axis(ax,'tight')

            %Label the axes
            label_shapespace_axes(ax,[],plot_info.stretch);

            %Set the tic marks
            set_tics_shapespace(ax,s);


            %Put an outline box around the plot
            nicebox(ax,'on')

            % 		%make hidden lines visible
            % 		hidden2('off',ax)

            %%%%
            %Make clicking on the thumbnail open it larger in a new window

            if ~plot_info.own_figure

                %build a plot_info structure for just the current plot
                plot_info_specific.axes = 'new';
                plot_info_specific.components = plot_info.components(i);
                plot_info_specific.category = 'beta';
                plot_info_specific.style = plot_info.style;
                plot_info_specific.CCFtype = plot_info.CCFtype;
                plot_info_specific.stretch = plot_info.stretch;

                %set the button down callback on the plot to be sys_draw with
                %the argument list for the current plot, and set the button
                %down callback for the mesh to the same
                set(plot_info.axes(i),'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});
                set(meshhandle,'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});

            else

                set(get(ax(1),'Parent'),'Name',[beta_list{function_number} ' Beta'])

                switch plot_info.style

                    case 'surface'

                        udata.plottype = 'CCF-surface';

                    case 'contour'

                        udata.plottype = 'CCF-contour';

                end

                %Mark this figure as betas
                udata = get(plot_info.figure(i),'UserData');
                udata.plottype = 'beta';
                set(plot_info.figure(i),'UserData',udata);

            end


        end
    else
        for i = 1:length(plot_info.axes)

            %call up the relevant axis
            ax =plot_info.axes(i);
            
            text(0.1,0.5,'Illustration of Beta not implemented for >2 shape variables','Parent',ax)
        end
    end

end