function plot_info = mfield_draw(s,p,plot_info,sys,shch,resolution)

% 6/23/17 - Copy of vfield_draw edited by Jacquelin to plot the metric
%   ellipse field instead of vector field. Checkbox locations replaced BVI
% 6/28 - Finished making edits. Currently, code manually removes the
%   beginning and end rows and columns depending on whether or not the
%   metric stretch has been selected. Should probably improve this in the
%   future. 

    %Get the configuration file, and extract the Colorpath
	configfile = 'sysplotter_config';
    configfile = fullfile(fileparts(mfilename('fullpath')),'..',configfile);
	load(configfile,'Colorset');
    
% Get the number of dimensions
n_dim = numel(s.grid.eval);

%Ensure that there are figure axes to plot into, and create new windows
%for those axes if necessary
plot_info = ensure_figure_axes(plot_info);


% Pull out the stretch name
stretchnames = {'stretch','surface'};
if plot_info.stretch
    stretchname = stretchnames{plot_info.stretch};
end


if n_dim==2
    
    %Vector field list
    vfield_list = {'X'};

    % Get the number of dimensions
    % 	n_dim = numel(s.grid.eval);


    %%%%%%%
    % Get the vector field and interpolate into the specified grid

    % Extract the plotting grid
    grid = s.grid.metric_display;


    % metric for evaluating ellipse field
    M = s.metricfield.metric_display.content.metric;
    %if ~plot_info.stretch
        tr = 2:(size(M{1},1)-1);
    %else
    %    tr = 3:(size(M{1},1)-2);
    %end
    % 
    % for i = 1:4
    %     
    %     Mtemp{i} = Mtemp{i}(tr,tr);
    %     
    % end
    %M = celltensorconvert(Mtemp);

    %%
    % for trimming  vectors that would be on the outside of the metric
    % stretch
    % Convert the ellipse field to the plotting grid specified in the gui
    resolution.metric = size(s.grid.metric_display{1});
    resolution.metric_range = [s.grid.metric_display{1}(1) s.grid.metric_display{1}(end) s.grid.metric_display{2}(1) s.grid.metric_display{2}(end)];
    [M,grid] = plotting_interp(M,grid,resolution,'metric');




    % Create a set of vector fields along coordinate directions
    E_coord = repmat({zeros(size(M{1}))},numel(grid));
    for i = 1:numel(grid)

        E_coord{i,i} = ones(size(M{1}));

    end

    %%%%%
    % If the shape coordinates should be transformed, make the conversion
    % (multiply the vectors by the inverse jacobian)

%     if plot_info.stretch
% 
%         M = celltensorconvert(M);
% 
%         % Calculate the jacobians at the plotting points
%         Jac = arrayfun(s.convert.stretch.jacobian,grid{:},'UniformOutput',false);
% 
%         % Use the jacobians to convert the metric
%         for i = 1:size(M,1)
%             for j = 1:size(M,2)
% 
%                 M{i,j} = (Jac{i,j}'\M{i,j})/Jac{i,j};
% 
%             end
%         end
% 
% 
%         % Use the jacobians to convert the coordinate vectors in accordance
%         % with the stretch transformation
% 
%         for i = 1:size(M{1},1)
% 
%             % Iterate over all vectors present
%             for j = 1:size(M{1},1)
% 
%                 % Extract all components of the relevant vector
%                 % 				tempEin = cellfun(@(x) x(j),E_coord(i,:));
%                 tempMin = [M{1}(i,j);M{2}(i,j)];
% 
%                 % Multiply by the Jacobian (because V_coord is a flow)
%                 tempMout = Jac{i,j}*tempMin;
% 
%                 % Replace vector components
% 
%                 M{1}(i,j) = tempMout(1);
%                 M{2}(i,j) = tempMout(2);
% 
%             end
% 
% 
%         end
% 
% 
%         % Convert the grid points to their new locations
%         [grid{:}] = s.convert.stretch.old_to_new_points(grid{:});
% 
%         M = celltensorconvert(M);
% 
%     end

    %%%
    %If there's a singularity, use arctan scaling on the magnitude of the
    %vector field
    if s.singularity

        M = arctan_scale_vector_fields(M);

    end


    % trim outer columns and rows

    for idx = 1:numel(M)
        M{idx} = M{idx}(tr,tr);
    end

    for idx = 1:numel(grid)
        grid{idx} = grid{idx}(tr,tr);
    end

    for i = 1:length(plot_info.axes)

        %call up the relevant axis
        ax = plot_info.axes(i);

        %get which vector field to use
        field_number = find(strcmp(plot_info.components{i}, vfield_list));

        %plot the vector field arrows

        % metricellipsefield(s.grid.metric_display{:},celltensorconvert(s.metricfield.metric_display.content.metric),'tissot',{'edgecolor','k’})

        %metricellipsefield(grid{:},celltensorconvert(M),'tissot-cross',{'edgecolor','k','parent',ax},{'color',Colorset.secondary,'parent',ax});
        metricellipsefield_convert(grid{:},celltensorconvert(M),'tissot-cross',s.convert,plot_info.stretch,{'linewidth',1,'edgecolor',Colorset.spot,'parent',ax},{'color',Colorset.secondary,'parent',ax});
        
        
        box(ax,'on');

        % Make edges and create a backing surfaceif coordinates have changed
        if plot_info.stretch

            edgeres = 30;

            oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
                s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
            oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
                linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];

            [x_edge,y_edge,z_edge] = s.convert.(stretchname).old_to_new_points(oldx_edge,oldy_edge);

            l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'ZData',z_edge,'Color','k','LineWidth',1); %#ok<NASGU>

            if plot_info.stretch == 2
                hold(ax,'on')
                [s_x,s_y,s_z] = s.convert.(stretchname).old_to_new_points(s.grid.eval{:});
                s_backing = surf('Parent',ax,'XData',s_x,'YData',s_y,'ZData',s_z,'FaceColor','w','EdgeColor','none');
                hold(ax,'off')
            end
            
        end


        if plot_info.stretch
            axis(ax,'equal');
            % 			axis(ax,[min(grid{1}(:)) max(grid{1}(:)) min(grid{2}(:)) max(grid{2}(:))]);
            axis(ax,[min(x_edge) max(x_edge) min(y_edge) max(y_edge)])
            if plot_info.stretch == 2
                view(ax,3)
            end
        else
            axis(ax,'equal');
        end
        %set the display range
        if ~plot_info.stretch
            axis(ax,s.grid_range);
        end

        %Label the axes (two-dimensional)
        label_shapespace_axes(ax,[],plot_info.stretch);

        %Set the tic marks
        set_tics_shapespace(ax,s);%,s.convert);

        %If there's a shape change involved, plot it
        if ~strcmp(shch,'null')
            if n_dim==2
                switch plot_info.stretch 
                    case {0,1} % No stretch or 2-d stretch
                        switch plot_info.style
                            case 'contour'
                                overlay_shape_change_2d(ax,p,plot_info.stretch,s.convert);
                            case 'surface'
                                overlay_shape_change_3d_surf(ax,p,zdata{function_number,:},plot_info.stretch,s.convert,true);
                        end
                        
                    case 2 % Surface-embedded stretch
                        
                        overlay_shape_change_metricsurf(ax,p,s.convert.surface.old_to_new_points,Colorset)
                        %overlay_shape_change_3d_surf(ax,p,grid_extra,plot_info.stretch,s.convert,false)
                        %overlay_shape_change_2d(ax,p,plot_info.stretch,s.convert);
                end
                        
            end
            if n_dim>2
                meshhandle.FaceAlpha=0.9;
                line('Parent',ax,'XData',p.phi_locus_full{i}.shape(:,1),'YData',p.phi_locus_full{i}.shape(:,2),'ZData',p.phi_locus_full{i}.shape(:,3),'Color',Colorset.spot,'LineWidth',6,'parent',ax);
            end
        end

        %%%%
        %Make clicking on the thumbnail open it larger in a new window

        if ~plot_info.own_figure

            %build a plot_info structure for just the current plot
            plot_info_specific.axes = 'new';
            plot_info_specific.components = plot_info.components(i);
            plot_info_specific.category = 'mfield';
            plot_info_specific.stretch = plot_info.stretch;

            %set the button down callback on the plot to be sys_draw with
            %the argument list for the current plot
            set(plot_info.axes(i),'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});

        else

            set(get(ax(1),'Parent'),'Name','Metric Ellipse Field')

            %Mark this figure as a metric ellipse field
            udata = get(plot_info.figure(i),'UserData');
            udata.plottype = 'mfield';
            set(plot_info.figure(i),'UserData',udata);


        end


    end
else
        for i = 1:length(plot_info.axes)

            %call up the relevant axis
            ax = plot_info.axes(i);
            text(0.1,0.5,'Illustration of metric field not implemented for >2 shape variables','Parent',ax)
            
        end
end

end