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
s.geometry.baseframe = 'com-mean';

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
% 
% directions_possible = ff2n(nlinks);
directions_key = ["SSS","SSR","SRS","SRR","RSS","RSR","RRS","RRR","AXES"];
%% Calculate the consistency of the entire 4D space

[consistent_system,consistency_count,consistent_A] = full_consistency_space_3link(s, adot_scale);

%% This script confirms slopes on the slope plot will be constant
% squeeze(consistent_A(1,1,1,1,:,:)) will get you an individual A
% j and k refer to a1 and a2
% I want to confirm that for each of the j and k combinations, A is equal
% for all of the l and m values (when in the same subsystem)
for j = 1:size(consistent_A,1)
    for k = 1:size(consistent_A,2)
        for subsystem_to_check = 1:8
            previous_A = zeros(3,2);
            mismatch_count = 0;
            for l = 1:size(consistent_A,3)
                for m = 1:size(consistent_A,4)
                    if consistent_system(j,k,l,m) == subsystem_to_check
                        A = squeeze(consistent_A(j,k,l,m,:,:));
                        if any(abs(A - previous_A) > 0.0001)
                            mismatch_count = mismatch_count + 1;
                        end
                        previous_A = A;
                    end
                end
            end
            if mismatch_count > 1
                [j k subsystem_to_check]
            end
        end
    end
end

%% An algorithm to get the consistent system at each angle of shapechange
% gradient(squeeze(consistent_system(5,5,:,:)))

% first get a circle of shapechange values
shape = [1; 1];
n = 100;
T = linspace(0,2*pi,n);
X = zeros(n); Y = zeros(n); Z = zeros(n);
C = zeros(n);
for i = 1:n
    theta = T(i);
    shapechange = [cos(theta); sin(theta)];
    X(i) = shapechange(1);
    Y(i) = shapechange(2);
    % and determine its subsystem and A value
    [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
    [subsystem, A] = determine_subsystem(s,shape,shapechange,J_full);
    x_slope = A(1,:) * shapechange;
    C(i) = subsystem;
    Z(i) = x_slope;
end

max_subsystem = max(C,[],'all');
min_subsystem = min(C(C~=0),[],'all');

surf(X,Y,Z,C)

colormap(jet(max_subsystem - min_subsystem + 1))
colorbar
colorbar('Ticks',min_subsystem:max_subsystem,...
         'TickLabels', directions_key(min_subsystem:max_subsystem))


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

center = 0; amplitude = 0.5;
gait = {@(t) center + amplitude * [cos(t); -sin(t)], ...
                 @(t) amplitude * [-sin(t); -cos(t)]};

gait_cycles = 3;
frames_per_cycle = 100;
[consistent_with_SSS, consistent_subsystem, trajectory, t_solution] = when_do_scales_change_things(s,gait, gait_cycles, frames_per_cycle);
nframes = gait_cycles * frames_per_cycle;

s.physics.drag_ratio
% Say which subsystems appear
% ovserved_subsystems = directions_key(unique(consistent_subsystem))
% Find out the percentage of time spent in each one
[C,ia,ic] = unique(consistent_subsystem);
a_counts = accumarray(ic,1);
percent_in_direction = [directions_key(C)' 100*a_counts/sum(a_counts)]

%%
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
