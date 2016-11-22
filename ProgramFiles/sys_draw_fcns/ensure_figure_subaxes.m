%For plots which go into subfigures, create figure window with axes for
%subplots if they aren't provided
function plot_info = ensure_figure_subaxes(plot_info)

	%If axes aren't included in the plot_info, or the axes are the string
	%'new', generate sufficient new windows and keep track of their
	%axes
	if ~isfield(plot_info,'axes') || isempty(plot_info.axes) || strcmp(plot_info.axes,'new')

		%Make sure there's an empty matrix in the axes field
		plot_info.figure = [];
		plot_info.axes = [];

		%Make a new figure window
		plot_info.figure(end+1,1) = figure;

		% 		%Set an empty background color and make the emptiness
		% 		%exportable
		% 		set(f,'InvertHardCopy','off','Color','none');
		
		% Create nicer figures for rows of two or three full-width figures
		if any(length(plot_info.components) == [2 3])
			
			switch length(plot_info.components)
				
				case 2
					
					axispos = [0.13 0.6 0.775 0.35;
						0.13 0.15 0.775 0.35];
					
				case 3
		
					axispos = [0.13 0.69 0.775 0.25;
						0.13 0.42 0.775 0.25;
						0.13 0.15 0.775 0.25];
					
			end
			
			%fill in the axes field
			for i = 1:length(plot_info.components)

				%create an axis in that window and store its handle
				plot_info.axes(end+1,1) = axes('Position',axispos(i,:),'Parent',plot_info.figure);

			end
			
		else

			%fill in the axes field
			for i = 1:length(plot_info.components)

				%create an axis in that window and store its handle
				plot_info.axes(end+1,1) = subplot(length(plot_info.components),1,i,'Parent',plot_info.figure);

			end
			
		end
		
		%mark as in its own figure window
		plot_info.own_figure = true;

	else

		%mark as being in the thumbnail window
		plot_info.own_figure = false;
	end

end