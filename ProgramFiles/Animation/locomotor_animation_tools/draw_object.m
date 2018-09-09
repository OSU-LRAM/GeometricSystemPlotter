function object = draw_object(object)
%draw an object, keeping track of its handles for edge and fill

%update location calculations
object = update_vertices(object);

% % First check whether object has moved since it was last drawn. If not,
% % return the object untouched.
% % Otherwise, store the current position to be used as "last position" in
% % future calls.
% if ~isempty(object.ref_pos.last.draw);
%     
%     if (object.ref_pos.now == object.ref_pos.last.draw);
% 
%         return;
% 
%     end
%     
% else
%     
%     object.ref_pos.last.draw = object.ref_pos.now;
%     
% end


% % Next, if there's already a handle for this object, update the 
% if ~isempty(object.graphics.handle) && ishandle(object.graphics.handle);
%         
%     delete(object.graphics.handle(find(object.graphics.handle ~= 0 & ...
%         ishandle(object.graphics.handle))));
%     
% end

%%%%
% Fill if specified

if ~isempty(object.graphics.fill)
    
	if ~isempty(object.graphics.handle) && all(ishandle(object.graphics.handle))
		
        if iscell(object.geom.vertices)
            for idx = 1:numel(object.geom.vertices)
                set(object.graphics.handle(1,idx),'XData',object.geom.vertices{idx}(:,1),...
                    'YData',object.geom.vertices{idx}(:,2),...
                    object.graphics.fill{:},'Parent',object.parent);
                
            end
        else

            set(object.graphics.handle(1,1),'XData',object.geom.vertices(:,1),...
                'YData',object.geom.vertices(:,2),...
                object.graphics.fill{:},'Parent',object.parent);
            
        end
		
	else
	
		%draw the polygon and get its handle
        
        if iscell(object.geom.vertices)
            for idx = 1:numel(object.geom.vertices)
                object.graphics.handle(1,idx) = ...
                    patch('XData',object.geom.vertices{idx}(:,1),'YData',object.geom.vertices{idx}(:,2),...
                        object.graphics.fill{:},'Parent',object.parent,object.graphics.edge{:});
            end
        else
        
            object.graphics.handle(1,1) = ...
                patch('XData',object.geom.vertices(:,1),'YData',object.geom.vertices(:,2),...
                    object.graphics.fill{:},'Parent',object.parent,object.graphics.edge{:});
        end
			
	end
    
end

% %%%%
% % Plot the edge if it is specified
% 
% if ~isempty(object.graphics.edge);
%     
% 	if ~isempty(object.graphics.handle) && all(ishandle(object.graphics.handle)) && numel(object.graphics.handle)>=2
% 		
% 		set(object.graphics.handle(1,2),'Xdata',object.geom.vertices(:,1),...
% 			'Ydata',object.geom.vertices(:,2),...
%             object.graphics.edge{:},'Parent',object.parent);
% 		
% 	else
% 		%draw the polygon edge and get its handle
% 
% 		object.graphics.handle(1,2) = ...
% 			line('Xdata',object.geom.vertices(:,1),'Ydata',object.geom.vertices(:,2),...
% 				object.graphics.edge{:},'Parent',object.parent);
% 	end
%     
% end

    
% %%%%
% % Plot the sense line if it is specified
% 
% if isfield(object.geom,'sense');
%     
%     %draw the polygon edge and get its handle
%         
%     object.graphics.handle(1,3) = ...
%         line('Xdata',object.geom.sense(:,1),'Ydata',object.geom.sense(:,2),...
%             object.graphics.sense{:});
%     
% end
% 
% if isfield(object.geom,'probe');
%     
%     %draw the three probe lines
%     
%     %draw the two closest rays in red, and the third in green
%     ray_colors = {'r','r','g'};
%     
%     for i = 1:3
%         object.graphics.handle(i,4) = ...
%             line('Xdata',[object.ref_pos.now(1) object.geom.probe(i,1)]...
%             ,'Ydata',[object.ref_pos.now(2) object.geom.probe(i,2)]...
%             ,'Color',ray_colors{i});
%     
%     
%     
% %             line( [xy(1),closest_points(i,1)], [xy(2),closest_points(i,2)]...
% %             ,'Color',ray_colors{i} );
%     end
% end

drawnow;