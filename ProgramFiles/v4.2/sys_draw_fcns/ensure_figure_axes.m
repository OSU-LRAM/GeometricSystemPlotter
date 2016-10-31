%For plots which are one plot per figure, create figures with axes if they
%don't already exist
function plot_info = ensure_figure_axes(plot_info)

    %If axes aren't included in the plot_info, or the axes are the string
    %'new', generate sufficient new windows and keep track of their
    %axes
    if ~isfield(plot_info,'axes') || isempty(plot_info.axes) || strcmp(plot_info.axes,'new')
        
        %Make sure there's an empty matrix in the axes field
        plot_info.figure = [];
		plot_info.axes = [];
        
        %fill in the axes field
        for i = 1:length(plot_info.components)
            
            %Make a new figure window
            plot_info.figure(end+1,1) = figure;
			
			%Set an empty background color and make the emptiness
			%exportable
			%MOVING THIS TO THE PLOT COMMANDS
			%set(f,'InvertHardCopy','off','Color','none');
            
            %create an axis in that window and store its handle
            plot_info.axes(end+1,1) = axes('parent',plot_info.figure(end,1));
            
        end
        
        %mark as in its own figure window
        plot_info.own_figure = true;
        
    else
        
        %mark as being in the thumbnail window
		if ~isfield(plot_info,'own_figure')
			plot_info.own_figure = false;
		end
    end

end