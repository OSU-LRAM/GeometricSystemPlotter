function overlay_shape_change_3d_surf(ax,p,zdata,convert)
%overlay shape changes onto the specified axis

	%plot all paths
	for i = 1:numel(p.phi_locus_full)

		if exist('convert','var') && ~isempty(convert)

		% Get the value by which to scale the height function
		ascale = arrayfun(@(x,y) 1/det(convert.jacobian(x,y)),p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2));

		% Apply the jacobian to the vectors
		zdata{i} = ascale.*zdata{i};

		[p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2)] ...
				= convert.old_to_new_points(p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2));



		end


		%draw the path itself
		patch('XData',p.phi_locus_full{i}.shape(:,1),'YData',p.phi_locus_full{i}.shape(:,2),'Zdata',zdata{i},'EdgeColor','k','FaceColor','none','LineWidth',6,'Parent',ax);

	end
    
end