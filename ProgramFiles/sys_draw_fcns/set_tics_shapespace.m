function set_tics_shapespace(ax,s)
%place the tic marks
    
    %set the tic fontsize
    set(ax,'FontSize',30,'FontName','Helvetica')
		
    set(ax,'XTick',s.tic_locs.x,'YTick',s.tic_locs.y)

end