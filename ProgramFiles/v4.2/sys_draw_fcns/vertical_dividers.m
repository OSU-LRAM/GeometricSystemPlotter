function vertical_dividers(ax,divide_x,y_min,y_max)
%separate and label the segments of multi-segment gaits

	if length(divide_x) >= 2
		%Draw dividing lines between segments
		for i = 1:length(divide_x)-1

			%the dividing line
			line('Parent',ax,'XData',[divide_x(i) divide_x(i)],'YData',[y_min y_max]...
				,'LineStyle','--'...
				,'Color','k'...
				,'LineWidth',2.5);

% 			%Location for the number
% 			divide_text_x = divide_x(i-1) + diff(divide_x(i-1:i))/2;
% 			divide_text_y = 3*(y_max + y_min)/10;

			%put in the number
			%             text(divide_text_x,divide_text_y,num2str(i-1)...
			%                 ,'FontName','Times','Interpreter','latex','FontSize',30 ...
			%                     ,'Color',[0 0 0],'HorizontalAlignment','center'...
			%                         ,'VerticalAlignment','middle','Parent',ax)


		end
	end


end