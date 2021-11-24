function plot_info = sys_draw(plot_structure,sys,shch,stretch2,progress,update,resolution,handles)

	%make sure plot data file is up to date
	if update
		sys_update(sys,shch,stretch2,progress,handles)
        set_plot_resolution_from_file(handles)
        
        %%%%%%%%%%%%%%
        % Update the desired vector and scalar plotting resolutions
        resolution.vector = str2num(get(handles.vectorresolution,'String'));
        resolution.scalar = str2num(get(handles.scalarresolution,'String'));
        resolution.vector_range = get(handles.vectorresolution,'UserData');
        resolution.scalar_range = get(handles.scalarresolution,'UserData');
    end

    % Get the setup configuration file
    configfile = './sysplotter_config';
    load(configfile,'datapath')

    %merge the system and shape change file names into one
    plot_data_file = [sys '__' shch];

    %load the system and path data
    load(fullfile(datapath, plot_data_file),'s','p')
    
    if s.conf_space == LieGroups.SE2
        %plot all plots called for
        for i = 1:length(plot_structure)

            %generate the handle of a specific plot command from the structure
            eval(['plot_command = @' plot_structure(i).category '_draw;']);

            %call that plot command
            plot_info(i,1) = plot_command(s,p,plot_structure(i),sys,shch,resolution);

        end
    else
        warning('Native plotting not yet supported for non-SE2 position spaces');
        plot_info = [];
    end

	%Show full progress bar
	waitbar2a(1,progress,'Finished plotting')
	
	% Refresh the run display
	refresh_runinfo_Callback([], [], handles);
	
end


% 

% 
% 

% 
% 

% 


% 
% 
% function [x_min,x_max,y_min,y_max] = ...
%     set_axis_limits(target_axis,xdata,ydata,xbuffer,ybuffer)
% %set the axis limits for the plot defined by xdata and ydata
% 
%     %set the x range
%     x_min = min(xdata(:)) - (xbuffer * range(xdata(:)));
%     x_max = max(xdata(:)) + (xbuffer * range(xdata(:)));
%     
%     %set the y range
%     y_min = min(ydata(:)) - (ybuffer * range(ydata(:)));
%     y_max = max(ydata(:)) + (ybuffer * range(ydata(:)));
%     
%     %if either in fact has no range, set a default range of +/- 1
%     if abs(x_max-x_min) < 10*eps
%         
%         x_center = x_min;
% 
%         x_min = x_center-1;
%         x_max = x_center+1;
% 
%     end
%     if abs(y_max-y_min) < 10*eps
%         
%         y_center = y_min;
% 
%         y_min = y_center-1;
%         y_max = y_center+1;
% 
%     end
%     
%     %apply the axis limits
%     axis(target_axis,[x_min,x_max,y_min,y_max]);
% 
% 
% end
% 
% 

%             
% 
% 

% 
% % function title_axis(target_axis,title_text)
% % %write and format the plot title
% % 
% %     %Formatting options
% %     format_list = {'FontName','Times','FontSize',24,'Interpreter','latex'};
% %     
% %     %write the title
% %     %title(target_axis,title_text,format_list{:});
% %     
% % end
% 
% 
% 
% 
% 
% 
% 

% 
% 
% 
% 
% 
% function xy = jacobian_inv_multiplier_helper(J,vx,vy)
% 
% 	% Multiply the inverse of the jacobian times the vector
% 	xy_temp = J\[vx;vy];
% 	
% 	% rearrange to be in the 3rd dimension direction
% 	xy(1,1,:) = xy_temp;
% 	
% end
