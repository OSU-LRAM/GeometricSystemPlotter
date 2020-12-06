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

direction_names = ["FF","BF","FB","BB"];

s = sysf_two_link_lowRe;
nlinks = length(s.geometry.linklengths);

s = ensure_connection_and_metric(s); % required for create_grids
s = create_grids(s)
a_grid = s.grid.eval{1};
adot_grid = s.grid.eval{1}*0.1;

%% test for consistency (count, combination of a and adot, how many sub-systems are consistent: have the link
% velocities according to that sub-system match up with the link velocities
% that sub-system is supposed to be the active sub-system.) We expect one
% sub-system to be consistent at each a-adot combination. Zero being
% consistent somewhere would be weird and concerning. More than one might
% happen in more complicated systems (but is observed to not happen here--
% unless a or adot are zero. When that happens, all consistent sub-systems
% result in the same body velocity so it doesn't matter, just pick any).
% TODO: Check if there's anything to document within
% identify_piecewise_system if there's time. I think it's all boring stuff
% so I'm skipping for now.
% Figures/results from this secion:
% Consistency Count: The Z axis displays how many sub-systems are consistent.
% There is exactly one consistent sub-system for each nonzero combination
% of a and adot. The a and adot axes are consistent with all four
% sub-systems. I think it works out that all four systems give the same
% velocities at these values.
% Consistency Summary: The Z axis is a key corresponding to a specific
% sub-system. Because (other than the a and adot axes) there is only one
% sub-system consistent at each point, the currently consistent one can be
% displayed. As explained in the colorbar, BF is dark blue, FB is cyan, and
% BB is yellow. FF doesn't show up, and BB wouldn't show up either, but
% because it's last and the code is simple, it's displayed on the axes.
% This plot shows that only sub-systems FB and BF are ever active, and it's
% on a simple, per-quadrant basis. The FB quadrants correspond to when the
% swimmer is "closing" (alpha is getting farther from zero. the shape of
% the swimmer is getting farther from straight.) and the BF quadrants
% correspond to when the swimmer is "opening". This conclusion found
% computationally here is confirmed with a symmetry proof (explained
% elsewhere), so we know it's true for all friction coefficients that
% maintain the assumptions in that proof.
% This result, that the choice of what sub-system to use is so simple and
% ends up using only two sub-systems instead of four, makes implementing
% and making sense of the full, scaled system so much easier. It shows that
% it was really worth it to start with such a simple and symmetric swimmer.
%%%%%%%%%%%%%%%%%%%

% find which system is consistent when
[system_map, count] = identify_piecewise_system(s, adot_grid);

% look at symmetry and missing/redundant areas
figure(3);
clf(3);
title("Consistency Count");
surface(adot_grid, a_grid, count);
xlabel("\alpha dot");
ylabel("\alpha");

figure(2);
clf(2);
title("Consistency Summary");
xlabel("\alpha dot");
ylabel("\alpha");
surface(adot_grid, a_grid, system_map);
colorbar('Ticks',[1, 2, 3, 4],...
         'TickLabels',direction_names)
    
%% demonstrate getting velocities and plotting as vectors in body and link
% origins
%%%%%%%%%%%%%%%%%%%

% get a body velocity
i = 20; j = 30;
bvel = apply_piecewise_system(s, system_map, 2, 0.1);

% get link velocities
lvel = zeros(nlinks, 3); % 3 is x y theta

% get Jfull(alpha)
[~, ~, J_full, ~, ~] = N_link_chain(s.geometry, a_grid(i));

for link = 1:nlinks
    lvel(link, :) = J_full{link} * [bvel; adot_grid(j)];
end

% plot
velocity_diagram(s, lvel, bvel, a_grid(i), adot_grid(j), 0, 0);

%% example of making some gaits
%%%%%%%%%%%%%%%%%%%

% define gait functions
centered_gait = generate_1D_gait(1, 0, 0);
offset_gait = generate_1D_gait(1, 1, 0);
offset_gait2 = generate_1D_gait(1, -1, 0);
tiny_gait = generate_1D_gait(0.5, 0, 0);
tiny_gait_bw = generate_1D_gait(0.5, 0, pi);
huge_gait = generate_1D_gait(3, 0, 0);
huge_gait_bw = generate_1D_gait(3, 0, pi);

gait = centered_gait;

%% find the displacement and trajectory the ODE solver way
%%%%%%%%%%%%%%%%%%%

% apply ODE solver
sol = asym_solve_gait(s, gait, system_map);
gait_displacement = deval(sol, 2*pi)
figure(3);
trajectory_plot(sol);

%% plot frames and turn it into a movie (takes a while)
%%%%%%%%%%%%%%%%%%%

F = animate_asymmetric_solution(s, sol, gait);

%% play the generated movie
%%%%%%%%%%%%%%%%%%%

% play movie
figure('visible','on') %forces animation to be visible in live script
movie(F);

%% save the movie (path might need adjusting)
%%%%%%%%%%%%%%%%%%%
% write movie
%v = VideoWriter('AsymmetricFriction/offset_gait.avi');
%open(v);
%writeVideo(v, F);
%close(v);

%% get A as a function of alpha (split into A.positive and A.negative
% depending on the sign of alphadot)
%%%%%%%%%%%%%%%%%%%

A = get_A_of_alpha(s, a_grid);

%% Make a plot for each of x, y, and theta components of A. Plot A+ and A-
% on the same axes.
%%%%%%%%%%%%%%%%%%%

% directional glyphs
figure()
tiledlayout(3,1)

% X plot
ax1 = nexttile;
plot(ax1,a_grid, A.positive(:,1), 'k>', a_grid, A.negative(:,1), 'r<');
ylabel(ax1,'A_x')
xlabel(ax1,'\alpha')

% Y plot
ax1 = nexttile;
plot(ax1,a_grid, A.positive(:,2), 'k>', a_grid, A.negative(:,2), 'r<');
ylabel(ax1,'A_y')
xlabel(ax1,'\alpha')

% Theta plot
ax2 = nexttile;
plot(ax2,a_grid, A.positive(:,3), 'k>', a_grid, A.negative(:,3), 'r<');
ylabel(ax2,'A_\theta')
xlabel(ax2, '\alpha')

% directional glyphs
figure()
tiledlayout(2,1)

% X plot
ax1 = nexttile;
plot(ax1,a_grid, A.difference(:,1),'k')
ylabel(ax1,'A_x')
xlabel(ax1, '\alpha')

% Theta plot
ax2 = nexttile;
plot(ax2,a_grid, A.difference(:,3),'k')
ylabel(ax2,'A_\theta')
xlabel(ax2, '\alpha')

%% Make slope plots in both styles. I think the "true" option with arrows on
% the positive direction might be nice.
%%%%%%%%%%%%%%%%%%%

A_sparse = get_A_of_alpha(s, a_grid(1:2:31));
slope_plot(a_grid(1:2:31),A_sparse,true)
slope_plot(a_grid(1:2:31),A_sparse,false)

%% compare a bunch of tiny gaits with different offsets: experimenting with
% lie brackets, essentially. Probably not something to put into the paper,
% but hey I did it.
%%%%%%%%%%%%%%%%%%%

amp = 0.1;
step = 0.1;
off_grid = -pi+amp : step : pi-amp;

displacements = zeros([length(off_grid),3]);
for i = 1:length(off_grid)
    off = off_grid(i);
    gait = generate_1D_gait(amp, off, 1);
    sol = asym_solve_gait(s, gait, system_map);
    displacements(i,:) = deval(sol, 2*pi);
end
% figure()
% tiledlayout(3,1)
% ax = nexttile;
% plot(ax, off_grid, displacements(:,1))
% ax = nexttile;
% plot(ax, off_grid, displacements(:,2))
% ax = nexttile;
% plot(ax, off_grid, displacements(:,3))

% the backwards version:
displacements_bw = zeros([length(off_grid),3]);
for i = 1:length(off_grid)
    off = off_grid(i);
    gait = generate_1D_gait(0.5, off, -1);
    sol = asym_solve_gait(s, gait, system_map);
    displacements_bw(i,:) = deval(sol, 2*pi);
end
% figure()
% tiledlayout(3,1)
% ax = nexttile;
% plot(ax, off_grid, displacements_bw(:,1))
% ax = nexttile;
% plot(ax, off_grid, displacements_bw(:,2))
% ax = nexttile;
% plot(ax, off_grid, displacements_bw(:,3))

% compare forwards and then backwards to backwards and then forwards
difference = zeros([length(off_grid),3]);
for i = 1:length(off_grid)
    fw = vec_to_mat_SE2(displacements(i, :));
    bw = vec_to_mat_SE2(displacements_bw(i, :));
    fw_then_bw = fw * bw;
    bw_then_fw = bw * fw;
    difference(i,:) = mat_to_vec_SE2(fw_then_bw) - mat_to_vec_SE2(bw_then_fw);
end
amp
figure()
tiledlayout(3,1)
ax = nexttile;
plot(ax, off_grid, difference(:,1))
ax = nexttile;
plot(ax, off_grid, difference(:,2))
ax = nexttile;
plot(ax, off_grid, difference(:,3))

amp
figure()
tiledlayout(3,1)
ax = nexttile;
plot(ax, off_grid, difference(:,1))
ax = nexttile;
plot(ax, off_grid, difference(:,2))
ax = nexttile;
plot(ax, off_grid, difference(:,3))


%% Looking at a spread of gaits, holding one variable constant, and seeing
% how the trajectory changes.
%%%%%%%%%%%%%%%%%%%

% testing out trajectory plots
figure()
hold on

amp = 0.5;
for off = -2:0.5:2

gait = generate_1D_gait(amp, off, 0);
sol = asym_solve_gait(s, gait, system_map);
trajectory_plot(sol);
end

hold off
axis equal

% comparing trajcetories of same gait, different phase
figure()
hold on

for phase = 0:pi/4:2*pi
    
    gait = generate_1D_gait(1.5, 0, phase);
    sol = asym_solve_gait(s, gait, system_map);
    trajectory_plot(sol)
    
end

hold off
axis equal

%% Experimenting with integrating Adif to get a prediction of the net
% displacement produced by a gait.
%%%%%%%%%%%%%%%%%%%

Ax_dif = A.difference(:,1);
Ay_dif = A.difference(:,2);
At_dif = A.difference(:,3);

point_a = a_grid(14)
point_b = a_grid(18)
da = a_grid(2) - a_grid(1)

dx = trapz(Ax_dif(14:18))*da
dy = trapz(Ay_dif(14:18))*da
dth = trapz(At_dif(14:18))*da

% compare to a gait that goes from point_a to point_b:
off = (point_a + point_b)/2;
amp = abs(point_a - point_b)/2;
gait = generate_1D_gait(amp,off,0);
sol = asym_solve_gait(s, gait, system_map);

displacement = deval(sol, 2*pi)'

g_cerc_vector = [dx dy dth];
g_circ = [0 -g_cerc_vector(3) g_cerc_vector(1); ...
          g_cerc_vector(3) 0 g_cerc_vector(2); ...
          0 0 0];
Adif_prediction = mat_to_vec_SE2(expm(g_circ))

error = Adif_prediction - displacement
proport_error = error ./ displacement

% test a bigger amp now that I think it's correct:
point_a = a_grid(11)
point_b = a_grid(21)
da = a_grid(2) - a_grid(1);

dx = trapz(Ax_dif(11:21))*da;
dy = trapz(Ay_dif(11:21))*da;
dth = trapz(At_dif(11:21))*da;

% compare to a gait that goes from point_a to point_b:
off = (point_a + point_b)/2;
amp = abs(point_a - point_b)/2;
gait = generate_1D_gait(amp,off,0);
sol = asym_solve_gait(s, gait, system_map);

displacement = deval(sol, 2*pi)'

g_cerc_vector = [dx dy dth];
g_circ = [0 -g_cerc_vector(3) g_cerc_vector(1); ...
          g_cerc_vector(3) 0 g_cerc_vector(2); ...
          0 0 0];
Adif_prediction = mat_to_vec_SE2(expm(g_circ))

error = Adif_prediction - displacement
proport_error = error ./ displacement

%% Using the function that packages this calculation into something less
% cumbersome
%%%%%%%%%%%%%%%%%%%

%Adif_gait_prediction(s, 1, 0, calc_density)
% this function is unfinished