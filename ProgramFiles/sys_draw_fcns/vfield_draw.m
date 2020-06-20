function plot_info = vfield_draw(s,p,plot_info,sys,shch,resolution)

    %Get the configuration file, and extract the Colorpath
	configfile = 'sysplotter_config';
    configfile = fullfile(fileparts(mfilename('fullpath')),'..',configfile);
	load(configfile,'Colorset');

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
    
    % Set all "small" values of V to zero, to avoid spurious arrows
    % appearing
    for idx = 1:size(V,1)
        
        normV = sqrt(V{idx,1}.^2 + V{idx,2}.^2);
        
        V{idx,1}(normV < 1e-10) = 0;
        V{idx,2}(normV < 1e-10) = 0;
        
    end
	
	%%
	% Convert the vector field to the plotting grid specified in the gui
	[V,grid] = plotting_interp(V,grid,resolution,'vector');
    
    % Create a set of vector fields along coordinate directions
    V_coord = repmat({zeros(size(V{1}))},numel(grid));
    for i = 1:numel(grid)

        V_coord{i,i} = ones(size(V{1}));

    end
	
	%%%%%
	% If the shape coordinates should be transformed, make the conversion
	% (multiply the vectors by the inverse jacobian)
	
	if plot_info.stretch && (n_dim == 2)
        
        if plot_info.stretch==1
            % Calculate the jacobians at the plotting points
            Jac = arrayfun(s.convert.stretch.jacobian,grid{:},'UniformOutput',false);
        end
        
        if plot_info.stretch==2
            % Calculate the jacobians at the plotting points
            Jac = arrayfun(s.convert.surface.jacobian,grid{:},'UniformOutput',false);
        end
		
		% Use the jacobians to convert the vectors
		V_converted = repmat({V{1,1}},size(V,1),size(Jac{1},1));
        
		% Iterate over all connection vector fields present
        for i = 1:size(V,1)
			
			% Iterate over all vectors present
            for j = 1:numel(V{i,1})
				
				% Extract all components of the relevant vector
				tempVin = cellfun(@(x) x(j),V(i,:));
				
				% Multiply by the inverse transpose Jacobian (because V is a
				% gradient, not a flow)
                
                % Expand Jacobian if mapping 2 to 3 dimensions. Expansion
                % is a vector normal to both columns, length does not
                % matter
                if size(Jac{j},1) == 3
                    J_exp = [Jac{j} cross(Jac{j}(:,1),Jac{j}(:,2))];
                    tempVin = [tempVin(:);0];
                else
                    J_exp = Jac{j};
                end
                
				tempVout = (J_exp')\tempVin(:);
				
				% Replace vector components
                for k = 1:numel(tempVout)
					
					V_converted{i,k}(j) = tempVout(k);
					
                end
				
            end
            
			
        end
        V = V_converted;


		% Use the jacobians to convert the coordinate vectors in accordance
		% with the stretch transformation
		    % Create a set of vector fields along coordinate directions
        V_coord = repmat({zeros(size(V{1}))},size(Jac{1},1));
        for i = 1:numel(grid)

            V_coord{i,i} = ones(size(V{1}));

        end
        V_coord_converted = V_coord;
		
		% Iterate over all coordinate vector fields present
        for i = 1:size(V_coord,1)
			
			% Iterate over all vectors present
            for j = 1:numel(V_coord{i,1})
				
				% Extract all components of the relevant vector
				tempVin = cellfun(@(x) x(j),V_coord(i,1:numel(grid)));
				
				% Multiply by the Jacobian (because V_coord is a flow)
				tempVout = Jac{j}*tempVin(:);
				
				% Replace vector components
				for k = 1:numel(tempVout)
					
					V_coord_converted{i,k}(j) = tempVout(k);
					
				end
				
            end
            
			
        end
        V_coord = V_coord_converted;
        
        
        % Rotate the coordinate vectors to get normals
        V_norm = repmat({zeros(size(V{1}))},numel(grid));
        if size(V_coord,2) == 2
            
            V_norm = {V_coord{2,2} -V_coord{2,1};
                -V_coord{1,2} V_coord{1,1}};
                
        elseif size(V_coord,2) == 3
            
            V_norm = {V_coord{3,3} V_coord{3,2} -V_coord{3,1}; % z coord field around y
                -V_coord{1,2} V_coord{1,1} -V_coord{1,3}; % x around z
                V_coord{2,1} V_coord{2,3} -V_coord{2,2}}; % y around x
            
        else
            warning('Rotations for >3 dimensions undefined in vfield_draw clipping of stretched vector field (add them if you need them) ')
            
        end
        
        % Convert the grid points to their new locations
		
        if plot_info.stretch==1
            [grid{:}] = s.convert.stretch.old_to_new_points(grid{:});
        end
        
        if plot_info.stretch==2
            [grid{:},grid_extra] = s.convert.surface.old_to_new_points(grid{:});
            grid = [grid;{grid_extra}];
        end
        
		
	
    
        %%%%
        % Trim any vectors that go outside the boundary
        if plot_info.stretch < 2
            % Take the dot product of the connection vector fields and the
            % normal vector fields
            dprods = repmat({zeros(size(V_norm{1}))},size(V,1),size(V_norm,1));
            for i = 1:size(dprods,1)

                for j = 1:size(dprods,2)

                    elementprods = cellfun(@(x,y) x.*y,V(i,:),V_norm(j,:),'UniformOutput',false);

                    for k = 1:numel(elementprods)

                        dprods{i,j} = dprods{i,j}+elementprods{k};

                    end

                end

            end

            % Create a cell array to hold the masking term
            edgemask = repmat({zeros(size(V{1}))},size(V));

            % Iterate along the rows of V (each row is one field)
            for idxA = 1:size(V,1)

                % Iterate along the elements of dotprods in the corresponding row
                % (each colum corresponds to the dot product of the current row of
                % V with the then nth coordinate field
                for idxB = 1:size(dprods,2)

                    % Identify the index sets that correspond to the first and last
                    % elements along this direction of the grid
                    indices_start = [repmat({':'},1,idxB-1), {1}, repmat({':'},1,size(V,2)-idxB)];
                    indices_end = [repmat({':'},1,idxB-1), {size(V{1},idxB)}, repmat({':'},1,size(V,2)-idxB)];

                    % Find all vectors that point out
                    V_test_start = dprods{idxA,idxB} < 0;
                    V_test_end = dprods{idxA,idxB} > 0;

                    % Take the start and end indices values from the test boolean
                    edgemask{idxA,idxB}(indices_start{:}) = V_test_start(indices_start{:});
                    edgemask{idxA,idxB}(indices_end{:}) = V_test_end(indices_end{:});

                end

                % Combine all the edgemasks for a field into a single mask

                edgemask_merged = zeros(size(edgemask{1}));
                for idxB = 1:size(V,2)

                    edgemask_merged = edgemask_merged | edgemask{idxA,idxB};

                end

                % Apply edgemask_merged to all the fields in this row of V 
                for idxB = 1:size(V,2)

                    V{idxA,idxB}(edgemask_merged) = 0;

                end       
            end
        end    
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
        field_number = find(strcmp(plot_info.components{i}, vfield_list));
        
		%plot the vector field arrows
        if n_dim == 2
            
            if numel(grid) == 2
                q = quiver(ax,grid{:},V{field_number,1},V{field_number,2},'k','LineWidth',2);
            else
                q = quiver3(ax,grid{:},V{field_number,1},V{field_number,2},V{field_number,3},'k','LineWidth',2);
            end
            
            
        else
            idxt=cell(1,n_dim-3);
            idxt(1,:)={1};
			q = quiver3(ax,grid{1}(:,:,:,idxt{:}),grid{2}(:,:,:,idxt{:}),grid{3}(:,:,:,idxt{:}),V{field_number,1}(:,:,:,idxt{:}),V{field_number,2}(:,:,:,idxt{:}),V{field_number,3}(:,:,:,idxt{:}),'k','LineWidth',2);
        end
			
        % Make edges if coordinates have changed
        if plot_info.stretch && (numel(s.grid.eval) == 2)

            edgeres = 30;

            oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
                s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
            oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
                linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];

            switch plot_info.stretch
                
                case 1 % 2-d stretch
            
                    [x_edge,y_edge] = s.convert.stretch.old_to_new_points(oldx_edge,oldy_edge);

                    l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Color','k','LineWidth',1);

                case 2 % 3-d surface embedding
            
                    [x_edge,y_edge,z_edge] = s.convert.surface.old_to_new_points(oldx_edge,oldy_edge);

                    l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Zdata',z_edge,'Color','k','LineWidth',1);
                    
                    hold(ax,'on')
                    [s_x,s_y,s_z] = s.convert.surface.old_to_new_points(s.grid.eval{:});
                    s_backing = surf('XData',s_x,'YData',s_y,'ZData',s_z,'Parent',ax,'FaceColor','w','EdgeColor','none');
                    hold(ax,'off')
            end

       end				
		
		if plot_info.stretch
			axis(ax,'equal');
			axis(ax,[min(grid{1}(:)) max(grid{1}(:)) min(grid{2}(:)) max(grid{2}(:))]);
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
        
        box(ax,'on')
        
        %Label the axes (two-dimensional)
        label_shapespace_axes(ax,[],plot_info.stretch);
		
        %Set the tic marks
        set_tics_shapespace(ax,s);
        
        %If there's a shape change involved, plot it
        if ~strcmp(shch,'null')
            if n_dim==2
                switch plot_info.stretch 
                    case {0,1} % No stretch or 2-d stretch
                        overlay_shape_change_2d(ax,p,plot_info.stretch,s.convert);
                        
                    case 2 % Surface-embedded stretch
                        
                        overlay_shape_change_metricsurf(ax,p,s.convert.surface.old_to_new_points,Colorset)
                        %overlay_shape_change_3d_surf(ax,p,grid_extra,plot_info.stretch,s.convert,false)
                        %overlay_shape_change_2d(ax,p,plot_info.stretch,s.convert);
                end
                        
            end
            if n_dim>2
                line('Parent',ax,'XData',p.phi_locus_full{i}.shape(:,1),'YData',p.phi_locus_full{i}.shape(:,2),'ZData',p.phi_locus_full{i}.shape(:,3),'Color',Colorset.spot,'LineWidth',6,'parent',ax);
            end
        end

        
        
        %%%%
        %Make clicking on the thumbnail open it larger in a new window
        
        if ~plot_info.own_figure

            %build a plot_info structure for just the current plot
            plot_info_specific.axes = 'new';
            plot_info_specific.components = plot_info.components(i);
            plot_info_specific.category = 'vfield';
			plot_info_specific.stretch = plot_info.stretch;
			
            %set the button down callback on the plot to be sys_draw with
            %the argument list for the current plot
            set(plot_info.axes(i),'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});
            
        else
            
            set(get(ax(1),'Parent'),'Name',[vfield_list{field_number} ' Vector Field'])
			
			%Mark this figure as a vector field
			udata = get(plot_info.figure(i),'UserData');
			udata.plottype = 'vfield';
			set(plot_info.figure(i),'UserData',udata);

        
        end
        
        
    end
    
end