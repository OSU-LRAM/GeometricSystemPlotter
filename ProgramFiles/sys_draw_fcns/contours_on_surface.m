function contours_on_surface(convert,h,grid)
    
    % Set the x, y, and z bounds as tight axes
    axis(h.Parent,...
        [min(grid{1}(:)) max(grid{1}(:)) min(grid{2}(:)) max(grid{2}(:)) min(grid{3}(:)) max(grid{3}(:))])
    
   

    
    % loop over all the contours, using undocumented
    % contour handle features
    for idx = 1:numel(h.EdgePrims)

        v = h.EdgePrims(idx).VertexData;
        [cx,cy,cz] = convert.surface.old_to_new_points(v(1,:),v(2,:));
        c = [cx;cy;cz];

        h.EdgePrims(idx).VertexData = single(c);

    end

    try
        for idx = 1:numel(h.EdgeLoopPrims)

            v = h.EdgeLoopPrims(idx).VertexData;
            [cx,cy,cz] = convert.surface.old_to_new_points(v(1,:),v(2,:));
            c = [cx;cy;cz];

            h.EdgeLoopPrims(idx).VertexData = single(c);

        end
    catch
    end
    
    for idx = 1:numel(h.FacePrims)
        
            v = h.FacePrims(idx).VertexData;
            [cx,cy,cz] = convert.surface.old_to_new_points(v(1,:),v(2,:));
            c = [cx;cy;cz];

            h.FacePrims(idx).VertexData = single(c);

        h.FacePrims(idx).ColorData = uint8([255 255 255 255]');
        
    end

    
end