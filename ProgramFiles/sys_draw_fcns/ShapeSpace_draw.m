function plot_info = ShapeSpace_draw(s,p,plot_info,sys,shch,resolution)

    % If the request is to display in optimized coordinates, append the change-of-body-frame
    if strcmp(plot_info.components{1},'opt')
        
        % Ensure that the baseframe is a cell array, so that the change to
        % optimized coordinates can be appended
        if ~iscell(s.geometry.baseframe)
            s.geometry.baseframe = {s.geometry.baseframe};
        end
        
        % Append the transformation to optimized coordinates
        s.geometry.baseframe = [s.geometry.baseframe {sys}];
        
        % Change the display name
        displayname = 'Minimum-perturbation coordinates';
    else
        displayname = s.geometry.baseframe;
    end
    
    plot_info = ensure_figure_axes(plot_info);
    

    illustrate_shapespace(s,plot_info.axes)
    
        if ~plot_info.own_figure

			%build a plot_info structure for just the current plot
			plot_info_specific.axes = 'new';
			plot_info_specific.components = plot_info.components;
			plot_info_specific.category = 'ShapeSpace';
			plot_info_specific.style = plot_info.style;
			plot_info_specific.CCFtype = plot_info.CCFtype;
			plot_info_specific.stretch = plot_info.stretch;

            %set the button down callback on the plot to be sys_draw with
			%the argument list for the current plot, and set the button
			%down callback for the mesh to the same
			set(plot_info.axes,'ButtonDownFcn',{@sys_draw_dummy_callback,plot_info_specific,sys,shch,plot_info.stretch_name});

		else

			set(get(plot_info.axes,'Parent'),'Name',['Baseframe: ' displayname])

			%Mark this figure as a constraint curvature function
			udata = get(plot_info.figure,'UserData');
			
			switch plot_info.style
				
				case 'surface'
					
					udata.plottype = 'CCF-surface';
					
				case 'contour'
					
					udata.plottype = 'CCF-contour';
					
			end
			
			set(plot_info.figure,'UserData',udata);
			
		end


end
