function label_shapespace_axes(ax,z_text,converted)
%put the x and y labels on the plots. optional argument 'z_text' sets a
%label for the z axis

    %Default text labels
    x_text = '$\alpha_1$';
    y_text = '$\alpha_2$';
	
	% converted text labels
	if exist('converted','var') && converted
		x_text = '$r_{p1}$';
		y_text = '$r_{p2}$';
		
	end
    
    %Formatting options
    format_list = {'FontName','Times','Interpreter','latex','FontSize',30};
    
    %apply the label
    xlabel(ax,x_text,format_list{:});
    ylabel(ax,y_text,format_list{:});
    
    %if there's a third dimension, label it
    if exist('z_text','var') && ~isempty(z_text)
        
        zlabel(ax,z_text,format_list{:});
        
    end
        
    
end