function F = animate_asymmetric_solution(system, solution, gait, write)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

display.aspect_ratio = 0.05;
display.sharpness = 0.1;
display.scales = true;
display.points = 100;

gait_period = 2*pi;
gait_cycles = 3;
frames_per_cycle = 30;
nframes = gait_cycles * frames_per_cycle;

t_solution = softspace(0,gait_period * gait_cycles, nframes);
g_solution = deval(solution, mod(t_solution, gait_period));
gait_displacement = vec_to_mat_SE2(deval(solution, gait_period));

% H264
% PNG-like better for high contrast, jpegs are photos
% record movie
if write
    wobj = VideoWriter('movtest1.avi','Motion JPEG AVI');
    wobj.FrameRate = 10;                  % frames per second (video speed)
    open(wobj);
end


f_w = 1200; f_h = 600;
f = figure();
f.Position(3:4) = [f_w f_h];
f.Color = 'white';
ax1 = subplot(1,2,1);
% axis off
ax2 = subplot(1,2,2);
% axis off
subplot(ax1)
axis(ax1,'equal')
axis(ax2,'equal')
axis(ax1,'off')
axis(ax2,'off')
% axis equal
% axis equal
ax1.Position = [0 0 0.5 1];
ax2.Position = [0.5 0 0.5 1];
%ax = gca;
%ax.NextPlot = 'replaceChildren';
ax1.XLim = [-0.5,1];
ax1.YLim = [-0.5,1]*2*f_h/f_w - 0.25;
ax2.XLim = [-0.5,1];
ax2.YLim = [-0.5,1]*2*f_h/f_w - 0.25;
% xlim([-0.6, 1]);
% ylim([-0.5, 0.5]);

%%%% new stuff
if true
    trajectory = zeros(size(g_solution));
    
    
end

F = struct('cdata',[],'colormap',[]); % reset length of struct
F(nframes) = struct('cdata',[],'colormap',[]);
for frame = 1:nframes
    t = t_solution(frame);
    gaits_finished = floor(t/gait_period);
    
    accumulated_displacement = gait_displacement ^ gaits_finished;
    g = accumulated_displacement * vec_to_mat_SE2(g_solution(:, frame));
    
    if true %%% new stuff
        trajectory(:, frame) = mat_to_vec_SE2(g);
    end
    
    [B,~,~,~,scales_range] = fat_chain(system.geometry, gait{1}(t), display);
    B = g * B;
    %B(2,:) = B(2,:); % Offset so not overlapping with velocity diagram
    cla
    hold on
    if true %%%% new stuff
%         plot(ax1,trajectory(1,1:frame), trajectory(2, 1:frame),'LineWidth',0.8)
%         plot(ax2,trajectory(1,1:frame), trajectory(2, 1:frame),'LineWidth',0.8)
        
            % get a body velocity. The relevant sub-system is automatically chosen.
        bvel = apply_piecewise_system(system, gait{1}(t), gait{2}(t));

        % get link velocities
        lvel = zeros(length(system.geometry.linklengths), 3); % the 3 is for: x y theta

        % get Jfull(alpha)
        [~, ~, J_full, ~, ~] = N_link_chain(system.geometry, gait{1}(t));

        for link = 1:length(system.geometry.linklengths)
            lvel(link, :) = J_full{link} * [bvel; gait{2}(t)];
        end
%         velocity_diagram(system, lvel, bvel, gait{1}(t), gait{2}(t), 0, 0, 1);
    end
    
    
%      body_frame_axes = g * [0.1 0 1; 0 0.1 1]';
%      quiver([1 1]*trajectory(1,frame), [1 1]*trajectory(2,frame), body_frame_axes(1,:), body_frame_axes(2,:), '-r');
    
    body_frame_axes = [0.1 0 1; 0 0.1 1] * Adjinv(g);
%     quiver(ax1,[1; 1]*trajectory(1,frame), [1; 1]*trajectory(2,frame), body_frame_axes(:,1), body_frame_axes(:,2), '-k');
%     quiver(ax2,[1; 1]*trajectory(1,frame), [1; 1]*trajectory(2,frame), body_frame_axes(:,1), body_frame_axes(:,2), '-k');
    
%     patch(ax1,B(1,:),B(2,:),[1 0 0])
%     patch(ax2,B(1,:),B(2,:),[1 0 0])
    for link = 1:length(system.geometry.linklengths)
        if lvel(link,1) < 0
            color = [0 0 0];
            lwidth = 1;
        else
            color = [1 1 1]*0.7;
            lwidth = 0.5;
        end
%         patch(ax1,B(1,scales_range(link,1):scales_range(link,2)),B(2,scales_range(link,1):scales_range(link,2)),[1 0 0],'EdgeColor',color,'LineWidth',lwidth)
%         patch(ax2,B(1,scales_range(link,1):scales_range(link,2)),B(2,scales_range(link,1):scales_range(link,2)),[1 0 0],'EdgeColor',color,'LineWidth',lwidth)
    end
    
    subplot(ax1)
    cla(ax1)
    hold(ax1,'on')
    plot(ax1,trajectory(1,1:frame), trajectory(2, 1:frame),'LineWidth',0.8)
    quiver(ax1,[1; 1]*trajectory(1,frame), [1; 1]*trajectory(2,frame), body_frame_axes(:,1), body_frame_axes(:,2), '-k');
    patch(ax1,B(1,:),B(2,:),[1 0 0])
    patch(ax1,B(1,scales_range(link,1):scales_range(link,2)),B(2,scales_range(link,1):scales_range(link,2)),[1 0 0],'EdgeColor',color,'LineWidth',lwidth)
    hold(ax1,'off')
    
    subplot(ax2)
    cla(ax2)
    hold(ax2,'on')
    plot(ax2,trajectory(1,1:frame), trajectory(2, 1:frame),'LineWidth',0.8)
    quiver(ax2,[1; 1]*trajectory(1,frame), [1; 1]*trajectory(2,frame), body_frame_axes(:,1), body_frame_axes(:,2), '-k');
    patch(ax2,B(1,:),B(2,:),[1 0 0])
    patch(ax2,B(1,scales_range(link,1):scales_range(link,2)),B(2,scales_range(link,1):scales_range(link,2)),[1 0 0],'EdgeColor',color,'LineWidth',lwidth)
    hold(ax2,'off')
    
    drawnow
    F(frame) = getframe;
    
    if write
        cdata = print('-RGBImage','-r200');     % save image with '-r200' resolution
        frame_img = im2frame(cdata);              % convert image to frame
        writeVideo(wobj,frame_img);           % save frame into video
    end
end

if write
    close(wobj); 
end

end

