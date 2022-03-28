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
	
    grid_orig = grid; % Save the original grid for contour-hacking
	if plot_info.stretch && (numel(s.grid.eval) == 2)
			
        switch plot_info.stretch
            
            case 1 % This is the 2-d stretch
                
                % Get the value by which to scale the constraint curvature function
                ascale = arrayfun(@(x,y) 1/det(s.convert.stretch.jacobian(x,y)),grid{:});

                % Apply the jacobian to the vectors
                H = cellfun(@(x) x.*ascale,H,'UniformOutput',false);

                % Convert the grid points to their new locations
                [grid{:}] = s.convert.stretch.old_to_new_points(grid{:});
                
            case 2 % This is embedding the surface in 3-d space
                
                % Get the value by which to scale the constraint curvature function
                J_grid = arrayfun(@(x,y) s.convert.surface.jacobian(x,y),grid{:},'UniformOutput',false);
                J_det = cellfun(@(J) norm(cross(J(:,1),J(:,2))),J_grid);

                % Apply the jacobian to the vectors
                H = cellfun(@(x) x./J_det,H,'UniformOutput',false);

                % Convert the grid points to their new locations
                [grid{:},grid_extra] = s.convert.surface.old_to_new_points(grid{:});
                
        end
               
		
	end
	
	
	


    %Make the plots
    for i = 1:length(plot_info.axes)
        
        %call up the relevant axis
        ax =plot_info.axes(i);
        
        %get which constraint curvature function to use
        function_number = strcmp(plot_info.components{i}, CCF_list);
        
      
        % Handle 2-d plots differently from 3-d plots
        if n_dim == 2

            % Handle metric surface plots differently from
            % unstretched or metric stretched plots
            if plot_info.stretch == 2
                
                % Specific plotting commands based on the plot style
                % selected
                switch plot_info.style

                    % "surface" plot with metric surface should use the
                    % surface from the metric, with color taken from the
                    % ccf
                    case 'surface'
                        
                        meshhandle = surf(ax,grid{:},grid_extra,H{function_number});
                        colormap(Colorset.colormap); 
                        shading(ax,'interp')
                        axis(ax,'equal')
                        axis(ax,'tight')
                        view(ax,3)
                        
                    % ideally, "contour" would allow us to overlay a
                    % contour plot of one function on the surface of
                    % another, but this isnt easy. code below uses
                    % undocumented matlab features
                    case 'contour'
                        
                        [~, meshhandle] = contour(ax,grid_orig{:},H{function_number},7,'linewidth',2,'Fill','on');
                        drawnow
                        contours_on_surface(s.convert,meshhandle,[grid; {grid_extra}])
                        addlistener(meshhandle, 'MarkedClean', @(src,evnt)contours_on_surface(s.convert,meshhandle,[grid; {grid_extra}]));
                       
                        colormap(Colorset.colormap_contour); 
                        %shading(ax,'interp')
                        axis(ax,'equal')
                        axis(ax,'tight')
                        view(ax,3)
                   otherwise
                        
                        error('plot_info.style does not contain a valid option')
                        
                end
            
            % This "else" matches unstretched or 2-d stretched plots
            else
                
                % Specific plotting commands based on the plot style
                % selected
                switch plot_info.style

                    % "surface" plot with metric surface should use the
                    % surface from the metric, with color taken from the
                    % ccf
                    case 'surface'
                        
                        meshhandle = surf(ax,grid{:},H{function_number});
                        colormap(Colorset.colormap); 
                        axis(ax,'tight')
                       %shading(ax,'interp')
                        view(ax,3)                         

                    case 'contour'
                        
                         [~, meshhandle] = contour(ax,grid{:},H{function_number},7,'linewidth',2);
                         axis(ax,'equal')
                        axis(ax,'tight')
                          
                        colormap(Colorset.colormap_contour); 
                   otherwise
                        
                        error('plot_info.style does not contain a valid option')
                        
                end
                
            end

        % This "else" catches the case where the shape space has more than
        % two dimensions. The code figures out the plane intersecting the
        % maximum-norm point on the ccf in the direction that most
        % interacts with it; it badly needs commenting
        else

            %switch plot_info.style

                %case 'surface'

%                     % create an empty cell array with as many entries as
%                     % there are dimensions
%                     y=cell(1,n_dim);
%                     
%                     % populate this cell array with zeros
%                     y(:)={0};
%                     
%                     % Build up a (squared) norm of the constraint curvature
%                     % (using for now the simple in-coordinates norm, could
%                     % introduce the metric later)
%                     Hnorm = zeros(size(H{1}));
%                     for idx_normbuild = 1:size(H,2)
%                         Hnorm = Hnorm+H{function_number,idx_normbuild}.^2;
%                     end
%                     
%                     % Get the maximum value, and the index of that value,
%                     % for the current CCF
%                     [~, Imax] = max(Hnorm(:));
%                     [Imax_x,Imax_y,Imax_z] = ind2sub(size(H{function_number}),Imax);
%                     
%                     % Populate y with the location of that point
%                     for idx_maxpoint = 1:n_dim
%                         y{idx_maxpoint} = grid{idx_maxpoint}(Imax);
%                     end
%                     
%                     
%                     
                    % create a cell array with two fewer elements than
                    % there are dimensions
                    idxt=cell(1,n_dim-2);
                    
                    % populate this cell array with ones
                    idxt(1,:)={1};
%                     
%                     % For every dimension in the problem, create a copy of
%                     % the grid? Not sure why this is named
%                     % "interpstatecurvature"
%                     for j=1:1:n_dim
%                         interpstatecurvature{j}=grid{j,1};
%                     end
% 
%                     % Get the value of the curvature at the point y
%                     for j=1:n_dim*(n_dim-1)/2
%                         curvature(:,j)=interpn(interpstatecurvature{:},H{function_number,j},y{:},'cubic');
%                     end
%                     
%                     % Create a skew-symmetric matrix describing the
%                     % orientation of the curvature form
%                     B=[0,curvature(1),curvature(2);
%                         -curvature(1),0,curvature(n_dim);
%                         -curvature(2),-curvature(n_dim),0];
% 
%                     % Get the eigenvectors and values of the curvature
%                     % orientation
%                     [V,D]=eig(B);
%                     
%                     % Sort the eigenvalues and eigenvectors
%                     [d,ind] = sort(diag(D));
%                     Ds = D(ind,ind);
%                     Vs = V(:,ind);
%                  
% 
%                     % Get the real component of the eigenvector with the
%                     % largest eigenvalue
%                     X=real(Vs(:,end));
%                     
%                     % Get a vector that is right-hand-positive orthogonal
%                     % to the preceding vector
%                     Y=(Vs(:,end)-X)/(sqrt(-1));
% 
%                     % Normalize the X and Y vectors (to make up for the
%                     % fact that they lost length when mapped to real values
%                     Xnorm=X/(norm(X));
%                     Ynorm=Y/(norm(Y));

%                     normal=cross(Xnorm,Ynorm);
%                     y=zeros(1,n_dim)
%                     projnorm=y*normal;

                    pointonplane=zeros(1,n_dim);%y'-(y'-projnorm*normal);
                    
%                     % Test code
%                     Xnorm = [1;0;0];
%                     Ynorm = [0; 0; 1];
%                     y = {0 0 0};

                    [maxpoint, maxplane, Imax_x, Imax_y, Imax_z,interpstatecurvature] = CCF_maxpoint(H(function_number,:),grid);
                    Xnorm = maxplane(:,1);
                    Ynorm = maxplane(:,2);

                    % Pull out a 2d x and y grid
                    Xtemp=grid{1,1}(:,:,idxt{:});
                    Ytemp=grid{2,1}(:,:,idxt{:});

                    %arsize=length(grid{1,1}(:,1));
                    idxt2=cell(1,n_dim-3);
                    idxt2(1,:)={0};
                    
                    % Loop over points on the x y grid
                    zgrid_at_max = Xtemp(Imax_x,Imax_y)*Xnorm(3)+Ytemp(Imax_x,Imax_y)*Ynorm(3);
                    for m=1:1:size(Xtemp,1)
                        for j=1:1:size(Xtemp,2)
                            
                                % Rotate the grid so that it lines up with
                                % the curvature two-form
                                xgrid(m,j)=Xtemp(m,j)*Xnorm(1)+Ytemp(m,j)*Ynorm(1);%+pointonplane(1);
                                ygrid(m,j)=Xtemp(m,j)*Xnorm(2)+Ytemp(m,j)*Ynorm(2);%+pointonplane(2);
                                zgrid(m,j)=Xtemp(m,j)*Xnorm(3)+Ytemp(m,j)*Ynorm(3);%+pointonplane(3);
                                
%                                 % Raise or lower the grid so that the
%                                 % maxpoint is on the surface
%                                 zshift = maxpoint{3}-zgrid_at_max;%(zgrid_at_max-y{3});
%                                 zgrid(m,j)=zgrid(m,j) + zshift;

                        end
                        
                    end
                    
                    pointdist = (xgrid-maxpoint{1}).^2 + (ygrid-maxpoint{2}).^2 + (zgrid-maxpoint{3}).^2;
                    [~,mindistI] = min(pointdist,[],'all','linear');
                    
                    planeshift = [maxpoint{1} - xgrid(mindistI); maxpoint{2} - ygrid(mindistI); maxpoint{3} - zgrid(mindistI)];
                    
                    xgrid = xgrid+planeshift(1);
                    ygrid = ygrid+planeshift(2);
                    zgrid = zgrid+planeshift(3);
                    
                    for m=1:1:size(Xtemp,1)
                        for j=1:1:size(Xtemp,2)
                                
                                for k=1:2
                                    curvaturetemp(:,k)=interpn(...
                                        interpstatecurvature{:},...
                                        H{function_number,k},...
                                        xgrid(m,j),ygrid(m,j),zgrid(m,j),...
                                        idxt2{:},...
                                        'spline');
                                end
                                curvaturetemp(:,3)=interpn(interpstatecurvature{:},H{function_number,n_dim},xgrid(m,j),ygrid(m,j),zgrid(m,j),idxt2{:},...
                                    'spline');
                                curvatureproj(m,j)=curvaturetemp(1)*(Xnorm(1)*Ynorm(2)-Ynorm(1)*Xnorm(2))+curvaturetemp(2)*(Xnorm(1)*Ynorm(3)-Ynorm(1)*Xnorm(3))+curvaturetemp(3)*(Xnorm(2)*Ynorm(3)-Ynorm(2)*Xnorm(3));
                        end
                    end

                    meshhandle=pcolor(xgrid,ygrid,-curvatureproj,'Parent',ax);
                    meshhandle.ZData=zgrid;
                    

                    hold on

                    colormap(Colorset.colormap); 
                    shading(ax,'interp')
                    view(ax,3)
                    
%                 case 'contour'
%                     
%                     y=cell(1,n_dim);
%                     y(:)={0};
%                     idxt=cell(1,n_dim-2);
%                     idxt(1,:)={1};
%                     for j=1:1:n_dim
%                         interpstatecurvature{j}=grid{j,1};
%                     end
% 
%                     for j=1:n_dim*(n_dim-1)/2
%                         curvature(:,j)=interpn(interpstatecurvature{:},H{function_number,j},y{:},'cubic');
%                     end
%                     
%                     B=[0,curvature(1),curvature(2);-curvature(1),0,curvature(n_dim);-curvature(2),-curvature(n_dim),0];
% 
%                     [V,D]=eig(B);
%                     [d,ind] = sort(diag(D));
%                     Ds = D(ind,ind);
%                     Vs = V(:,ind);
%                     
%                     
%                     X=real(Vs(:,end));
%                     Y=(Vs(:,end)-X)/(sqrt(-1));
%                     
%                     Xnorm=X/(norm(X));
%                     Ynorm=Y/(norm(Y));
% 
%                     Xtemp=grid{1,1}(:,:,idxt{:});
%                     Ytemp=grid{2,1}(:,:,idxt{:});
% 
%                     arsize=length(grid{1,1}(:,1));
%                     idxt2=cell(1,n_dim-3);
%                     idxt2(1,:)={ceil(arsize/2)};
%                     idxt3=cell(1,n_dim-3);
%                     idxt3(1,:)={0};
%                     for m=1:1:arsize
%                         for j=1:1:arsize
%                                 xgrid(m,j)=Xtemp(m,j)*Xnorm(1)+Ytemp(m,j)*Ynorm(1);
%                                 ygrid(m,j)=Xtemp(m,j)*Xnorm(2)+Ytemp(m,j)*Ynorm(2);
%                                 zgrid(m,j)=Xtemp(m,j)*Xnorm(3)+Ytemp(m,j)*Ynorm(3);
%                         end
%                     end
% 
% 
%                     for idx1=1:length(grid{1,1}(:,1))
%                         for idx2=1:length(grid{2,1}(:,1))
%                             for idx3=1:length(grid{3,1}(:,1))
%                                 for k=1:n_dim*(n_dim-1)/2
%                                     curvaturetemp(:,k)=interpn(interpstatecurvature{:},H{function_number,k},grid{2,1}(idx1,idx2,idx3,idxt2{:}),grid{1,1}(idx1,idx2,idx3,idxt2{:}),grid{3,1}(idx1,idx2,idx3,idxt2{:}),idxt3{:},'cubic');
%                                 end                                
%                                 curvatureproj(idx1,idx2,idx3)=curvaturetemp(1)*(Xnorm(1)*Ynorm(2)-Ynorm(1)*Xnorm(2))+curvaturetemp(2)*(Xnorm(1)*Ynorm(3)-Ynorm(1)*Xnorm(3))+curvaturetemp(n_dim)*(Xnorm(2)*Ynorm(3)-Ynorm(2)*Xnorm(3));
%                             end
%                         end
%                     end
%                                
% 
%                    meshhandle=contourslice(grid{2,1}(:,:,:,idxt2{:}),grid{1,1}(:,:,:,idxt2{:}),grid{3,1}(:,:,:,idxt2{:}),curvatureproj,xgrid,ygrid,zgrid,'parent',ax);
%                    view(ax,3)                    
% 
%                     hold on
% 
%                     colormap(Colorset.colormap_contour); 
%                     view(ax,3)
%             end


        end
                        

				
				
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
                if strcmp(plot_info.style,'surface')
                    meshhandle.FaceAlpha=0.9;
                end
                line('Parent',ax,'XData',p.phi_locus_full{i}.shape(:,1),'YData',p.phi_locus_full{i}.shape(:,2),'ZData',p.phi_locus_full{i}.shape(:,3),'Color',Colorset.spot,'LineWidth',6,'parent',ax);
            end
        end

        %tight axes
        axis(ax,'tight')


        %Put an outline box around the plot
        box(ax,'on')
        set(ax,'BoxStyle','back');

        % Make edges if coordinates have changed
        if plot_info.stretch && (numel(s.grid.eval) == 2) %&& (strcmp(plot_info.style,'contour'))

            edgeres = 30;

            oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
                s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
            oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
                linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];

            switch plot_info.stretch
                
                case 1 % 2-d stretch
                    
                    if ~strcmp(plot_info.style,'surface')
            
                        [x_edge,y_edge] = s.convert.stretch.old_to_new_points(oldx_edge,oldy_edge);

                        l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Color','k','LineWidth',1);
                        
                    end

                case 2 % 3-d surface embedding
            
                    [x_edge,y_edge,z_edge] = s.convert.surface.old_to_new_points(oldx_edge,oldy_edge);

                    l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Zdata',z_edge,'Color','k','LineWidth',1);
            end
            
        elseif numel(s.grid.eval) > 2
            
            edgeres = 30;
            
            oldy_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
                s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
            oldx_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
                linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];
            
            x_edge = oldx_edge * Xnorm(1) + oldy_edge * Ynorm(1) +planeshift(1);
            y_edge = oldx_edge * Xnorm(2) + oldy_edge * Ynorm(2) +planeshift(2);
            z_edge = (oldx_edge * Xnorm(3)) + (oldy_edge * Ynorm(3)) +planeshift(3);
            
            l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Zdata',z_edge,'Color',[.5 .5 .5],'LineWidth',1);
       end
				

%         %equal axes sized to match grid or new dimensions if
%         %stretched
% 
%         if plot_info.stretch && (numel(s.grid.eval) == 2)
%             axis(ax,'equal');
%             axis(ax,[min(grid{1}(:)) max(grid{1}(:)) min(grid{2}(:)) max(grid{2}(:))]);
%         else
%             axis(ax,'equal','tight');
%         end
				
				
		

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

		%set(parH,'Colormap',coloration);

		%center the color map around zero
            Clim = get(ax,'Clim'); %get the current color limits
            C_outer = max(abs(Clim)); %get the maximum distance from zero
            set(ax,'Clim',[-C_outer C_outer]); %set an inclusive range around zero



		%Label the axes
		label_shapespace_axes(ax,[],plot_info.stretch);

		%Set the tic marks         
        set_tics_shapespace(ax,s);

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
			set(plot_info.axes(i),'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch});
			set(meshhandle,'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch});

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