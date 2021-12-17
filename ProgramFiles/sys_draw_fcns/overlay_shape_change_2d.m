function overlay_shape_change_2d(ax,p,stretch,convert,isomap,s)
%overlay shape changes onto the specified axis

    %Get the configuration file, and extract the Colorpath
	configfile = './sysplotter_config';
	load(configfile,'Colorset');

    %plot all paths
	for i = 1:numel(p.phi_locus)


		if stretch
            if stretch==1
                [p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2)] ...
                    = convert.stretch.old_to_new_points(p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2));
            else
                [p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2)] ...
                    = convert.surface.old_to_new_points(p.phi_locus_full{i}.shape(:,1), p.phi_locus_full{i}.shape(:,2));
            end

		end


		%draw the path itself
        
        for idx = 1:size(p.phi_locus_full{i}.shape,2)
            pathpoints{idx} = p.phi_locus_full{i}.shape(:,idx);
        end
        
%         if stretch==2
%         
%             H_path = griddata(isomap.x_new,isomap.y_new,isomap.HF_isomap,pathpoints{1},pathpoints{2});
% 
%             line(pathpoints{:},'ZData',H_path,'Color',Colorset.spot,'LineWidth',5,'Parent',ax);
% 
%   %          %draw the direction arrows on the line
%   %          plot_dir_arrows(p.phi_locus_full{i}.shape(:,1),p.phi_locus_full{i}.shape(:,2),p.phi_arrows{i}{1},isomap,s,'Color',Colorset.spot,'LineWidth',4,'Parent',ax);
% 
%         else
            
            line(pathpoints{:},'Color',Colorset.spot,'LineWidth',5,'Parent',ax);

            %draw the direction arrows on the line
            plot_dir_arrows(p.phi_locus_full{i}.shape(:,1),p.phi_locus_full{i}.shape(:,2),p.phi_arrows{i}{1},'Color',Colorset.spot,'LineWidth',4,'Parent',ax);

%         end

		%draw the start/end markers
		if ~isempty(p.phi_locus_full{i}.marker.shape)
			line(p.phi_locus_full{i}.marker.shape(:,1),p.phi_locus_full{i}.marker.shape(:,2),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',15,'Parent',ax,'LineStyle','none');
		end

		%If there are multiple marker points, use them
		if (max(size(p.phi_locus)) > 1) || (size(p.phi_locus_full{i}.marker.shape,1) > 1)

			for k = 1:size(p.phi_locus_full{i}.marker.shape,1) %draw all the markers for the path

				%Decide how to draw the current marker

				%If multiple paths in one config file
				if (max(size(p.phi_locus)) > 1)

					%If multiple markers for the current path
					if (size(p.phi_locus_full{i}.marker.shape,1) > 1)

						%make the marker for the first entry identify
						%the path number
						if k == 1

							marker_text = num2str(i);

							%the rest identify the section number
						else

							marker_text = num2str(k-1);

						end

					else

						%If only one section, identify the section number
						marker_text = num2str(i);

					end

				else

					%marker identifies section number
					marker_text = num2str(k);

				end





				%background_disk
				line(p.phi_locus_full{i}.marker.shape(k,1), p.phi_locus_full{i}.marker.shape(k,2),'Marker','o'...
					,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]...
					,'MarkerSize',30,'Parent',ax)

				%number
				text(p.phi_locus_full{i}.marker.shape(k,1), p.phi_locus_full{i}.marker.shape(k,2), marker_text...
					,'FontName','Helvetica','FontSize',30 ...
					,'Color',[1 0 0],'HorizontalAlignment','center'...
					,'VerticalAlignment','middle','Parent',ax)

			end

		end



	end
    
end