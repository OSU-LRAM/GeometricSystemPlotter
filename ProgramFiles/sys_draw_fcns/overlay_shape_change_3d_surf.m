function overlay_shape_change_3d_surf(ax,p,zdata,stretch,convert,scale_when_stretched)
%overlay shape changes onto the specified axis

    stretchnames = {'stretch'};

	%plot all paths
	for i = 1:numel(p.phi_locus_full)

		if stretch

		% Get the value by which to scale the z data
        if scale_when_stretched
            ascale = arrayfun(@(x,y) 1/det(convert.(stretchnames{stretch}).jacobian(x,y)),p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2));
        else
            ascale = 1;
        end
        
		% Apply the scale factor
		zdata{i} = ascale.*zdata{i};

		[p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2)] ...
				= convert.(stretchnames{stretch}).old_to_new_points(p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2));



		end


		%draw the path itself
		patch('XData',p.phi_locus_full{i}.shape(:,1),'YData',p.phi_locus_full{i}.shape(:,2),'Zdata',zdata{i},'EdgeColor','k','FaceColor','none','LineWidth',6,'Parent',ax);

	end
    
end