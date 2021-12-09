%% definitions
% Direction names (FF, BF...) refer to the direction of link velocity for
% each link. These "sub-systems" have constant drag matrices as long as you
% stay within one, so the old tools for scale-free systems apply. There are
% other things to call this that might be more
% intuitive, like "Smooth-Smooth" and "Rough-Smooth" which takes away the
% ambiguity of which thing is "forward" (link velocity or longitudinal
% friction?). s stores the information about the system, although not all
% of the fields are implemented. For example, I'm not worrying about the metric
% yet. Changing n_dim to 1 broke a lot of the existing code in the standard
% sysplotter software, which is why I wrote my own script to do the same
% sort of stuff. In s.geometry, you can see that the lengths are both 1. In
% s.physics, the drag properties are identical for each link. I use this
% link symmetry in the proof explaining why FB and BF are the only
% "sub-systems" that can happen, but I suspect it would still stand if the
% magnitude of drag numbers or lengths differed between links. It just
% simplifies things a lot.
% The rest of this code block is setting up the grids for the shape
% (referred to as "a" sometimes, but it corresponds to alpha) and shape
% velocity (referred to as "adot") to test for "consistency" (that for each
% combination of a and adot, there is exactly one sub-system where the link
% velocities according to that sub-system match up with the link velocities
% that sub-system is supposed to be the active sub-system.)
%%%%%%%%%%%%%%%%%%%

s = sysf_three_link_lowRe;
nlinks = length(s.geometry.linklengths);
%s.geometry.baseframe = 'com-mean';

% Define Drags
s.physics.drag_ratio = 2;
s.physics.drag_coefficient = 1;
s.physics.drag_bw_ratio = 4;
% FOR NULL CASE
% s.physics.drag_bw_ratio = 1;

s = ensure_connection_and_metric(s); % required for create_grids
s = create_grids(s);
a_grid = s.grid.eval{1};
adot_scale = 0.1;
adot_grid = s.grid.eval{1}*adot_scale;

directions_possible = ff2n(nlinks);
directions_key = ["SSS","SSR","SRS","SRR","RSS","RSR","RRS","RRR","AXES"];
%%

consistent_system = zeros(length(a_grid), length(a_grid), length(adot_grid), length(adot_grid));
consistency_count = zeros(length(a_grid), length(a_grid), length(adot_grid), length(adot_grid));

% loop through shapes
for j = 1:size(s.grid.eval{1},1)
    for k = 1:size(s.grid.eval{1},2)
        shape = [s.grid.eval{1}(j,k);s.grid.eval{2}(j,k)];
        [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
        % loop through every combination of Rough/Smooth links
        for i = 1:size(directions_possible,1)
            links_backwards = directions_possible(i,:);
            friction_directions = ones(1,nlinks) - 2 * links_backwards;
            [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry,s.physics, shape, links_backwards);
            % loop through shapechanges
            for l = 1:size(s.grid.eval{1},1)
                for m = 1:size(s.grid.eval{1},2)
                    shapechange = adot_scale*[s.grid.eval{1}(l,m);s.grid.eval{2}(l,m)];
                    body_vel = -A * shapechange;
                    
                    all_links_agree = 1;
                    for link = 1:length(s.geometry.linklengths)
                        link_body_vel = J_full{link} * [body_vel; shapechange];
                        dx_direction = sign(link_body_vel(1));
                        if(dx_direction == 0) % non-if way to do this?
                            agreement = 1; % if the link isn't moving, just say it agrees either way
                        else
                            agreement = friction_directions(link) * dx_direction;
                        end
                        all_links_agree = all_links_agree && (1 + agreement);
                    end
                    if all_links_agree
                        consistent_system(j, k, l, m) = i;
                        consistency_count(j, k, l, m) = consistency_count(j, k, l, m) + 1;
                    end
                    if consistency_count(j, k, l, m) == 8
                        % special code 9 for when all systems are
                        % consistent (axes)
                        consistent_system(j, k, l, m) = 9;
                    end
                end
            end
    end
    end
end

%%
% Trying to visualize the consistent_system four-dimensional array:


% Plotting slices with a constant shape:
shape_indices = [1,10];
shape = s.grid.eval{1}(shape_indices);
figure()
surface(adot_scale*s.grid.eval{1},adot_scale*s.grid.eval{2},reshape(consistent_system(shape_indices(1),shape_indices(2),:,:),size(s.grid.eval{1})))
colormap(colorcube(9))
colorbar
colorbar('Ticks',1:9,...
         'TickLabels',directions_key)
xlabel("d\alpha_1 / dt")
ylabel("d\alpha_2 / dt")
title("\alpha = ("+string(shape(1))+", "+string(shape(2))+")")
%
shape_indices = [17,29];
shape = s.grid.eval{1}(shape_indices);
figure()
surface(adot_scale*s.grid.eval{1},adot_scale*s.grid.eval{2},reshape(consistent_system(shape_indices(1),shape_indices(2),:,:),size(s.grid.eval{1})))
colormap(colorcube(9))
colorbar
colorbar('Ticks',1:9,...
         'TickLabels',directions_key)
xlabel("d\alpha_1 / dt")
ylabel("d\alpha_2 / dt")
title("\alpha = ("+string(shape(1))+", "+string(shape(2))+")")

% Plotting slices with a constant shapechange:
shapechange_indices = [1,10];
shapechange = adot_scale * s.grid.eval{1}(shapechange_indices);
figure()
surface(s.grid.eval{1},s.grid.eval{2},consistent_system(:,:,shapechange_indices(1),shapechange_indices(2)))
colormap(colorcube(9))
colorbar
colorbar('Ticks',1:9,...
         'TickLabels',directions_key)
xlabel("\alpha_1")
ylabel("\alpha_2")
title("d\alpha / dt = ("+string(shapechange(1))+", "+string(shapechange(2))+")")

% Plotting slices with a constant shapechange:
shapechange_indices = [20,20];
shapechange = adot_scale * s.grid.eval{1}(shapechange_indices);
figure()
surface(s.grid.eval{1},s.grid.eval{2},consistent_system(:,:,shapechange_indices(1),shapechange_indices(2)))
colormap(colorcube(9))
colorbar
colorbar('Ticks',1:9,...
         'TickLabels',directions_key)
xlabel("\alpha_1")
ylabel("\alpha_2")
title("d\alpha / dt = ("+string(shapechange(1))+", "+string(shapechange(2))+")")

%% Trying to do an animation
%
center = 0; amplitude = 1;
gait = {@(t) center + amplitude * [cos(t); sin(t)], ...
                 @(t) amplitude * [-sin(t); cos(t)]};
sol = asym_solve_gait(s, gait);

figure()
trajectory_plot(sol,3)

%% SSS on a standard gait: when is it FFF?
%

s = sysf_three_link_lowRe;
nlinks = length(s.geometry.linklengths);
s.geometry.baseframe = 'com-mean';

% Define Drags
s.physics.drag_ratio = 2; % vary this later
s.physics.drag_coefficient = 1;
s.physics.drag_bw_ratio = 1; % makes it always SSS

s = ensure_connection_and_metric(s); % required for create_grids
s = create_grids(s);
% a_grid = s.grid.eval{1};
% adot_scale = 0.1;
% adot_grid = s.grid.eval{1}*adot_scale;

directions_possible = ff2n(nlinks);
directions_key = ["SSS","SSR","SRS","SRR","RSS","RSR","RRS","RRR","AXES"];

center = 1; amplitude = 0.5;
gait = {@(t) center + amplitude * [cos(t); -sin(t)], ...
                 @(t) amplitude * [-sin(t); -cos(t)]};
sol = asym_solve_gait(s, gait);

% when is it going forward all links

t = 0; % time
shape = gait{1}(t);
directions_possible = ff2n(nlinks);
[~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
i = 1; % SSS
links_backwards = directions_possible(i,:);
friction_directions = ones(1,nlinks) - 2 * links_backwards;
% [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry,s.physics, shape% , links_backwards);
% body_vel = A * shapechange;


% trajectory_plot(sol,3)

gait_cycles = 3;
gait_period = 2*pi;
frames_per_cycle = 100;
nframes = gait_cycles * frames_per_cycle;

t_solution = linspace(0, gait_period*gait_cycles, nframes);
g_solution = deval(sol, mod(t_solution, gait_period));

gait_displacement = vec_to_mat_SE2(deval(sol, gait_period))


trajectory = zeros(size(g_solution));
consistency_on_path = zeros(1,nframes);
consistent_subsystem = zeros(1,nframes);

for frame = 1:nframes
    t = t_solution(frame);
    gaits_finished = floor(t/gait_period);
    
    accumulated_displacement = gait_displacement ^ gaits_finished;
    g = accumulated_displacement * vec_to_mat_SE2(g_solution(:, frame));
    
    
    trajectory(:, frame) = mat_to_vec_SE2(g);
    
    % now color code with the current consistency
    shape = gait{1}(t);
    [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
    shapechange = gait{2}(t);
    link_directions = [0 0 0];
    [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry,s.physics, shape);
    body_vel = -A * shapechange;
    for link = 1:length(s.geometry.linklengths)
        link_body_vel = J_full{link} * [body_vel; shapechange];
        dx_direction = sign(link_body_vel(1));
        link_directions(link) = dx_direction;
    end
    all_links_agree = all(~link_directions | ~(link_directions - friction_directions));
    if all(link_directions)
        % get the numbered subsystem currently being observed according to
        % link motion
        consistent_subsystem(frame) = bin2dec(num2str( (-link_directions+1)/2 )) + 1;
    else
        % if any link direction is currently zero, use a dummy value
        consistent_subsystem(frame) = 9;
    end
    
    
%     % redone above
%     for link = 1:length(s.geometry.linklengths)
%         link_body_vel = J_full{link} * [g_solution(:, frame); shapechange];
%         dx_direction = sign(link_body_vel(1));
%         if(dx_direction == 0) % non-if way to do this?
%             agreement = 1;
%         else
%             agreement = friction_directions(link) * dx_direction;
%         end_
%         all_links_agree = all_links_agree && (1 + agreement);
%     end
    
    consistency_on_path(frame) = all_links_agree;
    
end

s.physics.drag_ratio
% Say which subsystems appear
% ovserved_subsystems = directions_key(unique(consistent_subsystem))
% Find out the percentage of time spent in each one
[C,ia,ic] = unique(consistent_subsystem);
a_counts = accumarray(ic,1);
percent_in_direction = [directions_key(C)' 100*a_counts/sum(a_counts)]

max_subsystem = max(consistent_subsystem);
min_subsystem = min(consistent_subsystem);

% hold on
% plot(trajectory(1,:), trajectory(2, :));

% scatter(trajectory(1,:), trajectory(2, :), 100, 1 + consistency_on_path, 's', 'filled')
% colormap(jet(2))
% colorbar
% colorbar('Ticks',1:2,...
%          'TickLabels',["Not SSS", "SSS"])

% in world space
figure()
scatter(trajectory(1,:), trajectory(2, :), 100, consistent_subsystem, 's', 'filled')
colormap(jet(max_subsystem - min_subsystem + 1))
colorbar
colorbar('Ticks',min_subsystem:max_subsystem,...
         'TickLabels', directions_key(min_subsystem:max_subsystem))
xlabel("world frame x")
ylabel("world frame y")
title("drag ratio = " + s.physics.drag_ratio)
axis equal
     
% in shape space
figure()
frame_shapes = gait{1}(t_solution(1:nframes));
scatter(frame_shapes(1,:), frame_shapes(2,:), 100, consistent_subsystem, 's', 'filled')
colormap(jet(max_subsystem - min_subsystem + 1))
colorbar
colorbar('Ticks',min_subsystem:max_subsystem,...
         'TickLabels', directions_key(min_subsystem:max_subsystem))
xlabel("\alpha_1")
ylabel("\alpha_2")
title("drag ratio = " + s.physics.drag_ratio)
xlim([-pi pi])
ylim([-pi pi])

     
% hold off
% 
% title("Drag Ratio = " + s.physics.drag_ratio);
% 
axis square

% cd = single(consistency_on_path + 2);
% drawnow
% set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)

%%
% Attempt to compare com mean and centered baseframes. Are their link
% velocities the same as expected?
% c_m = sysf_three_link_lowRe;
% c_m.geometry.baseframe = 'com-mean';
% 
% % Define Drags
% c_m.physics.drag_ratio = 2; % vary this later
% c_m.physics.drag_coefficient = 1;
% c_m.physics.drag_bw_ratio = 1; % makes it always SSS
% 
% c_l = sysf_three_link_lowRe;
% %s.geometry.baseframe = 'com-mean';
% 
% % Define Drags
% c_l.physics.drag_ratio = 2; % vary this later
% c_l.physics.drag_coefficient = 1;
% c_l.physics.drag_bw_ratio = 1; % makes it always SSS
% 
% for frame = 1:nframes
%     t = t_solution(frame);
%     gaits_finished = floor(t/gait_period);
%     
%     % now color code with the current consistency
%     shape = gait{1}(t);
%     shapechange = gait{2}(t);
%     link_directions = [0 0 0];
%     [A_c_m, ~, ~, ~, ~] = LowRE_connection_discrete(c_m.geometry,c_m.physics, shape);
%     body_vel = -A * shapechange;
%     for link = 1:length(c_m.geometry.linklengths)
%         link_body_vel = J_full{link} * [body_vel; shapechange];
%     end
%     
% end

% FIXED I had failed to update J_Full!!!!