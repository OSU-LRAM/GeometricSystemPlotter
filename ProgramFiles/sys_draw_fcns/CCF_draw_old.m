function plot_info = CCF_draw_old(s,p,plot_info,sys,shch,resolution)
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
	
	if plot_info.stretch && (numel(s.grid.eval) == 2)
		
        if plot_info.stretch==1
            % Get the value by which to scale the constraint curvature function
            ascale = arrayfun(@(x,y) 1/det(s.convert.stretch.jacobian(x,y)),grid{:});

            % Apply the jacobian to the vectors
            H = cellfun(@(x) x.*ascale,H,'UniformOutput',false);

            % Convert the grid points to their new locations
            [grid{:}] = s.convert.stretch.old_to_new_points(grid{:});
        end
		
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
                if n_dim==2
				%Plot the constraint curvature function
                    meshhandle = surf(ax,grid{:},H{function_number});
                end
                
                 if n_dim>2
%                     curvinterest=H{function_number,:};
%                     for j=1:length(curvinterest(:,1,1))
%                         for k=1:length(curvinterest(1,:,1))
%                             for l=1:length(curvinterest(1,1,:))
%                                 totalcurvvalue(j,k,l)=sqrt(H{function_number,1}(j,k,l)^2+H{function_number,2}(j,k,l)^2+H{function_number,3}(j,k,l)^2);
%                             end
%                         end
%                     end
% 
%                     [max1,ind1]=max(totalcurvvalue);
%                     [max2,ind2]=max(max1);
%                     [max3,ind3]=max(max2);
% 
%                     totalcurvvalue(ind1(1,ind2(1,1,ind3),ind3),ind2(1,1,ind3),ind3);
% 
%                     y=[grid{1}(ind1(1,ind2(1,1,ind3),ind3),ind2(1,1,ind3),ind3),grid{2}(ind1(1,ind2(1,1,ind3),ind3),ind2(1,1,ind3),ind3),grid{3}(ind1(1,ind2(1,1,ind3),ind3),ind2(1,1,ind3),ind3)];
%                     
                    y=cell(1,n_dim);
                    y(:)={0};
                    idxt=cell(1,n_dim-2);
                    idxt(1,:)={1};
                    for j=1:1:n_dim
                        interpstatecurvature{j}=grid{j,1};
                    end

                    for j=1:n_dim*(n_dim-1)/2
                        curvature(:,j)=interpn(interpstatecurvature{:},H{function_number,j},y{:},'cubic');
                    end
                    
                    B=[0,curvature(1),curvature(2);-curvature(1),0,curvature(n_dim);-curvature(2),-curvature(n_dim),0];

                    [V,D]=eig(B);
                    [d,ind] = sort(diag(D));
                    Ds = D(ind,ind);
                    Vs = V(:,ind);
                    
                    
                    X=real(Vs(:,end));
                    Y=(Vs(:,end)-X)/(sqrt(-1));
                    
                    Xnorm=X/(norm(X));
                    Ynorm=Y/(norm(Y));

%                     normal=cross(Xnorm,Ynorm);
%                     y=zeros(1,n_dim)
%                     projnorm=y*normal;
                    
                    pointonplane=zeros(1,n_dim);%y'-(y'-projnorm*normal);
                    
                    Xtemp=grid{1,1}(:,:,idxt{:});
                    Ytemp=grid{2,1}(:,:,idxt{:});

                    arsize=length(grid{1,1}(:,1));
                    idxt2=cell(1,n_dim-3);
                    idxt2(1,:)={0};
                    for m=1:1:arsize
                        for j=1:1:arsize
                                xgrid(m,j)=Xtemp(m,j)*Xnorm(1)+Ytemp(m,j)*Ynorm(1)+pointonplane(1);
                                ygrid(m,j)=Xtemp(m,j)*Xnorm(2)+Ytemp(m,j)*Ynorm(2)+pointonplane(2);
                                zgrid(m,j)=Xtemp(m,j)*Xnorm(3)+Ytemp(m,j)*Ynorm(3)+pointonplane(3);
                                for k=1:2
                                    curvaturetemp(:,k)=interpn(interpstatecurvature{:},H{function_number,k},xgrid(m,j),ygrid(m,j),zgrid(m,j),idxt2{:},'spline');
                                end
                                curvaturetemp(:,3)=interpn(interpstatecurvature{:},H{function_number,n_dim},xgrid(m,j),ygrid(m,j),zgrid(m,j),idxt2{:},'spline');
                                curvatureproj(m,j)=curvaturetemp(1)*(Xnorm(1)*Ynorm(2)-Ynorm(1)*Xnorm(2))+curvaturetemp(2)*(Xnorm(1)*Ynorm(3)-Ynorm(1)*Xnorm(3))+curvaturetemp(3)*(Xnorm(2)*Ynorm(3)-Ynorm(2)*Xnorm(3));
                        end
                    end
                    
                    meshhandle=pcolor(xgrid,ygrid,-curvatureproj,'Parent',ax);
                    meshhandle.ZData=zgrid;
                    
                    hold on
                    
                    colormap(Colorset.colormap); 
                    shading(ax,'interp')
                    view(ax,3)
                    
                    
                end
				
				
				%If there's a shape change involved, plot it
				if ~strcmp(shch,'null')
                    if n_dim==2
					overlay_shape_change_3d_surf(ax,p,zdata{function_number,:},plot_info.stretch,s.convert,true);
                    end
                    if n_dim>2
                        meshhandle.FaceAlpha=0.9;
                        line('Parent',ax,'XData',p.phi_locus_full{i}.shape(:,1),'YData',p.phi_locus_full{i}.shape(:,2),'ZData',p.phi_locus_full{i}.shape(:,3),'Color',Colorset.spot,'LineWidth',6,'parent',ax);
                    end
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
				if n_dim==2
                    %Plot the constraint curvature function
                    [junk, meshhandle] = contour(ax,grid{:},H{function_number},7,'linewidth',2);				
                end
                
                if n_dim>2
%                     
%                     curvinterest=H{function_number,:};
% 
%                     for j=1:length(curvinterest(:,1,1))
%                         for k=1:length(curvinterest(1,:,1))
%                             for l=1:length(curvinterest(1,1,:))
%                                 totalcurvvalue(j,k,l)=sqrt(H{function_number,1}(j,k,l)^2+H{function_number,2}(j,k,l)^2+H{function_number,3}(j,k,l)^2);
%                             end
%                         end
%                     end
% 
%                     [max1,ind1]=max(totalcurvvalue);
%                     [max2,ind2]=max(max1);
%                     [max3,ind3]=max(max2);
% 
%                     totalcurvvalue(ind1(1,ind2(1,1,ind3),ind3),ind2(1,1,ind3),ind3)
% 
%                     y=[grid{1}(ind1(1,ind2(1,1,ind3),ind3),ind2(1,1,ind3),ind3),grid{2}(ind1(1,ind2(1,1,ind3),ind3),ind2(1,1,ind3),ind3),grid{3}(ind1(1,ind2(1,1,ind3),ind3),ind2(1,1,ind3),ind3)];
%                     
                     y=cell(1,n_dim);
                    y(:)={0};
                    idxt=cell(1,n_dim-2);
                    idxt(1,:)={1};
                    for j=1:1:n_dim
                        interpstatecurvature{j}=grid{j,1};
                    end

                    for j=1:n_dim*(n_dim-1)/2
                        curvature(:,j)=interpn(interpstatecurvature{:},H{function_number,j},y{:},'cubic');
                    end
                    
                    B=[0,curvature(1),curvature(2);-curvature(1),0,curvature(n_dim);-curvature(2),-curvature(n_dim),0];

                    [V,D]=eig(B);
                    [d,ind] = sort(diag(D));
                    Ds = D(ind,ind);
                    Vs = V(:,ind);
                    
                    
                    X=real(Vs(:,end));
                    Y=(Vs(:,end)-X)/(sqrt(-1));
                    
                    Xnorm=X/(norm(X));
                    Ynorm=Y/(norm(Y));

                    Xtemp=grid{1,1}(:,:,idxt{:});
                    Ytemp=grid{2,1}(:,:,idxt{:});

                    arsize=length(grid{1,1}(:,1));
                    idxt2=cell(1,n_dim-3);
                    idxt2(1,:)={ceil(arsize/2)};
                    idxt3=cell(1,n_dim-3);
                    idxt3(1,:)={0};
                    for m=1:1:arsize
                        for j=1:1:arsize
                                xgrid(m,j)=Xtemp(m,j)*Xnorm(1)+Ytemp(m,j)*Ynorm(1);
                                ygrid(m,j)=Xtemp(m,j)*Xnorm(2)+Ytemp(m,j)*Ynorm(2);
                                zgrid(m,j)=Xtemp(m,j)*Xnorm(3)+Ytemp(m,j)*Ynorm(3);
                        end
                    end


                    for idx1=1:length(grid{1,1}(:,1))
                        for idx2=1:length(grid{2,1}(:,1))
                            for idx3=1:length(grid{3,1}(:,1))
                                for k=1:n_dim*(n_dim-1)/2
                                    curvaturetemp(:,k)=interpn(interpstatecurvature{:},H{function_number,k},grid{2,1}(idx1,idx2,idx3,idxt2{:}),grid{1,1}(idx1,idx2,idx3,idxt2{:}),grid{3,1}(idx1,idx2,idx3,idxt2{:}),idxt3{:},'cubic');
                                end                                
                                curvatureproj(idx1,idx2,idx3)=curvaturetemp(1)*(Xnorm(1)*Ynorm(2)-Ynorm(1)*Xnorm(2))+curvaturetemp(2)*(Xnorm(1)*Ynorm(3)-Ynorm(1)*Xnorm(3))+curvaturetemp(n_dim)*(Xnorm(2)*Ynorm(3)-Ynorm(2)*Xnorm(3));
                            end
                        end
                    end
                               

                   meshhandle=contourslice(grid{2,1}(:,:,:,idxt2{:}),grid{1,1}(:,:,:,idxt2{:}),grid{3,1}(:,:,:,idxt2{:}),curvatureproj,xgrid,ygrid,zgrid,'parent',ax);
                   view(ax,3)

                    
                end
                                 
                    
				%If there's a shape change involved, plot it
				if ~strcmp(shch,'null')
                    if n_dim==2
					overlay_shape_change_2d(ax,p,plot_info.stretch,s.convert);
                    end
                    if n_dim>2
                        line('Parent',ax,'XData',p.phi_locus_full{i}.shape(:,1),'YData',p.phi_locus_full{i}.shape(:,2),'ZData',p.phi_locus_full{i}.shape(:,3),'Color',Colorset.spot,'LineWidth',6,'parent',ax);
                    end
    
				end
				
				
				% Make edges if coordinates have changed
				if plot_info.stretch && (numel(s.grid.eval) == 2)
					
					edgeres = 30;
					
					oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
						s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
					oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
						linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];

					[x_edge,y_edge] = s.convert.stretch.old_to_new_points(oldx_edge,oldy_edge);
					
					l_edge = line('Parent',ax,'Xdata',x_edge,'YData',y_edge,'Color','k','LineWidth',1);
					
				end
				
				%Put an outline box around the plot
				box(ax,'on')

				%equal axes sized to match grid or new dimensions if
				%stretched
				
				if plot_info.stretch && (numel(s.grid.eval) == 2)
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

                

                [x_new, y_new] = s.convert.surface.old_to_new_points(grid{1},grid{2});
                
%                 C = (s.convert.EI.C);
%                  C = reshape(smooth(C(:)),[],12);
                
                HF_isomap = griddata(s.convert.surface.EI.A,s.convert.surface.EI.B,s.convert.surface.EI.C,x_new, y_new,'cubic');
                
                isomap.x_new = x_new;
                isomap.y_new = y_new;
                isomap.HF_isomap = HF_isomap;
				
				%Plot the constraint curvature function
% 				[junk, meshhandle] = contour(ax,grid{:},H{function_number},7,'linewidth',2);
                meshhandle = pcolor(ax,grid{1},grid{2},H{function_number});
                
                meshhandle.ZData = HF_isomap;
                meshhandle.XData = x_new;
                meshhandle.YData = y_new;
                
                if plot_info.stretch==2
                    [grid{:}] = s.convert.surface.old_to_new_points(grid{:});
                end
                
                if plot_info.stretch==1
                    [grid{:}] = s.convert.stretch.old_to_new_points(grid{:});
                end
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

					if plot_info.stretch==2
                        [x_edge,y_edge] = s.convert.surface.old_to_new_points(oldx_edge,oldy_edge);
                    end
                    
                    if plot_info.stretch==1
                        [x_edge,y_edge] = s.convert.stretch.old_to_new_points(oldx_edge,oldy_edge);
                    end
                    
                    
                    HF_isomap_edge = griddata(s.convert.surface.EI.A,s.convert.surface.EI.B,s.convert.surface.EI.C,x_edge,y_edge);
					
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
				coloration = Colorset.colormap;
                
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