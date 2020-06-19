function overlay_shape_change_metricsurf(ax,p,conversion)


	%plot all paths
    for i = 1:numel(p.phi_locus_full)
        
        [x,y,z] = conversion(p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2));

        %draw the path itself
		patch('XData',x,'YData',y,'Zdata',z,'EdgeColor','k','FaceColor','none','LineWidth',6,'Parent',ax);
        
    end
    
end