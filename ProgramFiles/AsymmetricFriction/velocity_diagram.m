function velocity_diagram(s, link_velocity, body_velocity, a, adot, draw_axes, draw_blur, draw_velocity)
    display.aspect_ratio = 0.05;
    display.sharpness = 0.1;
    % convert link velocities into world coordinates
    [h, ~, ~, ~, ~] = N_link_chain(s.geometry, a);
    nlinks = length(h.lengths);
    l_vel_w = zeros(size(link_velocity));
    for link = 1:size(link_velocity,1)
        l_vel_w(link,:) = TeLg(h.pos(link,:))*link_velocity(link,:)';
    end
    if draw_axes
        x_axes = zeros(nlinks, 3);
        x_axes(:,1) = 1;
        for link = 1:size(x_axes,1)
            x_axes(link,:) = TeLg(h.pos(link,:))*x_axes(link,:)';
        end
        y_axes = zeros(nlinks, 3);
        y_axes(:,2) = 1;
        for link = 1:size(y_axes,1)
            y_axes(link,:) = TeLg(h.pos(link,:))*y_axes(link,:)';
        end
    end
    % now plot
    figure()
    axis equal
    % draw link outlines
    if draw_blur
        blurLength = 1;
        for t = 1:-0.2:0
             B = vec_to_mat_SE2(-body_velocity*blurLength*t) ...
               * fat_chain(s.geometry, a - adot*blurLength*t, display);
%              B = fat_chain(s.geometry, a - adot*blurLength*t, display);
%             B(:,1:102) = vec_to_mat_SE2(-l_vel_w(1,:)*blurLength*t) * B(:,1:102);
%             B(:,103:204) = vec_to_mat_SE2(-l_vel_w(2,:)*blurLength*t) * B(:,103:204);
            color = [1 1 1] -exp(-3*t*[1 1 1]);
            patch(B(1,:),B(2,:),[1 0 0],'EdgeColor',color )
            patch(B(1,:),B(2,:),[1 0 0],'EdgeColor',color )
        end
    else
        B = fat_chain(s.geometry, a, display);
        patch(B(1,:),B(2,:),[1 0 0])
    end
    hold on;
    % draw velocity arrows
    if draw_axes
        quiver(h.pos(:,1), h.pos(:,2), x_axes(:,1), x_axes(:,2), '--k');
        quiver(h.pos(:,1), h.pos(:,2), y_axes(:,1), y_axes(:,2), '--k');
    end
    if draw_velocity
        quiver(0, 0, body_velocity(1), body_velocity(2), 0.3/norm(body_velocity(1:2)),'b'); % body velocity
        quiver(h.pos(:,1), h.pos(:,2), l_vel_w(:,1), l_vel_w(:,2),'r'); % link velocities
    end
    hold off;
end

