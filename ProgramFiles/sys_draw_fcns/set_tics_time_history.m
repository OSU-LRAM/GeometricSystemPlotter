function set_tics_time_history(ax,s,timemarks)
%place the tic marks
    
    %set the tic fontsize
    set(ax,'FontSize',20,'FontName','Helvetica')
	
	% remove time marks if not called for
	if ~timemarks
		set(ax,'XTickLabel',[])
	end
	
end
