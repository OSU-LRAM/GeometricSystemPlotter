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

s.geometry.baseframe = 'com-mean';

% FOR NULL CASE
%s.physics.drag_bw_ratio = 1;

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
    
%% Demonstrate getting velocities and plotting as vectors in body and link
% origins
%%%%%%%%%%%%%%%%%%%

% Pick an a and adot on the grid. I can't remember if this is required for
% any of the code-- it might not be.
% i = 20; j = 30;
% a = a_grid(i)
% adot = adot_grid(j)
a = 2;
adot = 0.5;

% get a body velocity. The relevant sub-system is automatically chosen.
bvel = apply_piecewise_system(s, a, adot);

% get link velocities
lvel = zeros(nlinks, 3); % the 3 is for: x y theta

% get Jfull(alpha)
[~, ~, J_full, ~, ~] = N_link_chain(s.geometry, a);

for link = 1:nlinks
    lvel(link, :) = J_full{link} * [bvel; adot];
end

% plot
a = 0; adot = 0;
velocity_diagram(s, lvel, bvel, a, adot, 0, 0, 1);
velocity_diagram(s, lvel, bvel, a, adot, 0, 1, 1);
velocity_diagram(s, lvel, bvel, a, adot, 1, 0, 1);

%% Some example gaits.
% Input: generate_1D_gait(amplitude, center, phase). The parameters
% describe a sine wave.
% Output: two anonymous functions. gait{1} is a(t) and gait{2} is adot(t).
%%%%%%%%%%%%%%%%%%%

% define gait functions
centered_gait = generate_1D_gait(1, 0, 0);
offset_gait = generate_1D_gait(1, 1, 0);
offset_gait2 = generate_1D_gait(1, -1, 0);
tiny_gait = generate_1D_gait(0.5, 0, 0);
tiny_gait_bw = generate_1D_gait(0.5, 0, pi);
huge_gait = generate_1D_gait(3, 0, 0);
huge_gait_bw = generate_1D_gait(3, 0, pi);

% Replace centered_gait here with any of the above gaits or your own, and
% it will be the focus of the following few sections of the script.
gait = centered_gait;

gait = generate_1D_gait(1.5, 0, pi/2);

%% Find the displacement and trajectory using an ODE solver.
% 
%%%%%%%%%%%%%%%%%%%

% Wrapper that does all the details of using ode45. Assumes a gait's period
% is 2pi, which generate_1D_gait makes sure of.
sol = asym_solve_gait(s, gait);

% The gait's displacement is the value of the solution after one period has
% elapsed.
gait_displacement = deval(sol, 2*pi)

% figure(3);
% The gait's trajectory can be plotted with this function. It is hardcoded
% to plot three cycles.
% trajectory_plot(sol,3);

%% Plot frames and turn it into a movie.
% Uses softspace and does three cycles. The plots are saved to F, which can
% be played back or saved in the following two code blocks.
%%%%%%%%%%%%%%%%%%%
save_movie = false

F = animate_asymmetric_solution(s, sol, gait,save_movie);

%% play the generated movie
%%%%%%%%%%%%%%%%%%%

% play movie
%figure('visible','on') %forces animation to be visible in live script
movie(F,5,5);

%% save the movie (path might need adjusting)
%%%%%%%%%%%%%%%%%%%
% write movie
%v = VideoWriter('AsymmetricFriction/offset_gait.avi');
%open(v);
%writeVideo(v, F);
%close(v);

%% Get A as a function of alpha (split into A.positive and A.negative
% depending on the sign of alphadot). For a normal, no-scales system, A(alpha) is
% independant of adot. For this system, the combination of FB and BF
% sub-systems, A(alpha) depends on only the sign of adot. A.positive (or,
% A+) and A.negative (or, A-) correspond to the A that is active when adot
% is positive or negative, respectively. There is not a one-to-one
% correspondance between A+, A- and the A's corresponding to FB, BF. Going
% back to the consistency diagram that showed where FB and BF are valid in
% the a-adot space, you can see that:
% A+ uses [the A from the sub-system] FB for positve alpha and BF for
% negative alpha.
% A- uses BF for positive alpha and FB for negative alpha.
% The A's for sub-systems FB and BF are continuous with respect to alpha,
% but because of this hop between the two at alpha=0, A+ and A- are able to
% have a discontinuity in derivative at that point. This happens for the x
% components, as seen in the plot in the next code block. The fact
% that there isn't a discontinuity in the A's themselves is because the
% values for A's of FB and BF are identical at that point.
%%%%%%%%%%%%%%%%%%%

A = get_A_of_alpha(s, a_grid);

%% Make a plot for each of x, y, and theta components of A. Plot A+ and A-
% on the same axes.
% The first plot shows the components of A broken down into Ax, Ay, and
% Atheta. On each of these subplots, the black line with triangles pointing
% in the direction of increasing alpha is A+, and the red with triangles
% pointing towards decreasing alpha is A-. As the swimmer increases its
% joint angle alpha, its resulting body velocity scales according to the
% value of A+ at its current alpha value. As it decreases its joint angle
% alpha, the same thing happens but according to the current value of A-.
% Because we are concerned with closed gaits and there is only one degree
% of freedom, the swimmer must go back over any alpha value it has crossed
% again, but from the opposite directon. We can then try to approximate a
% single connection function "Adif" that is A+ - A-. Instead of
% "integrating" (not the exact thing going on if you're being fancy: slope
% plots are more) along the gait as we would do to get net displacement
% from the first figure, we integrate from the minimum alpha of a gait to
% the maximum alpha of a gait. As long as there's no extra backtracking.
% For example, a gait that was just the first gait x2 would go twice as
% far. There's a simpler way to explain this. The second plot shows Adif
% for x and theta. y is not shown because A+ and A- are equal everywhere.
% The only source of y displacement is from lie bracket trickery.
%%%%%%%%%%%%%%%%%%%

% directional glyphs
figure()
title("A- and A+");
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

figure()
title("Adif");
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
% The more mathematically appropriate way to view the connection function
% is as a field of slopes, one at each combination of shape variables (in
% a single-jointed snake, that's one at each value of alpha). (This has to
% do with it being co-vectors instead of vectors-- right? So scaling is
% different.) The value of A_x (or A_y, A_theta) is not a height on some
% axis, but the slope of the function relating alpha to x when alpha is at
% that point. For one shape variable, you can get away with viewing it as
% integration(?) but that perspective doesn't generalize correctly to
% higher dimensions of shape-space. A vector field might be a more helpful
% metaphor?
% The first way of plotting, with the multicolored X's, shows A+ and A-
% slope plots overlaid. The A value at the center of these line
% segments (where they cross in the x and theta subplots) determines the
% slope of the line segment. The width of the line segments, the change in
% alpha, is constant for each line segment. This way, the magnitude of the
% slope is clearly visible from the height reached by the line segment.
% However the A+ is only every used for increasing alpha and A- for
% decreasing alpha, so the part of the red line in the negative direction
% is never "used" (It's just a slope so it isn't actually being used, but
% if you think of it metaphorically as a ramp, it's a one-directional ramp
% from the center of the segment.) and vise-versa for A-. To make a more
% intuitive plot, I removed the "unused" portions to create the second
% plot, with the V's. Think of this one as a ramp that's "bent" in the
% middle. A system without scales would have one straight line for all of
% its ramps. This system has scales, which "bend" its one straight line.
% For example, on the Ax sub-plot, starting at any nonzero alpha and
% increasing alpha will push you up the red A+ ramp: A+x is positive for
% all nonzero alpha, so the slope is positive. This means the x portion of
% the body velocity will be increasing during this motion. Now, decreasing
% alpha puts us on the black A- ramp. The slope for A-x is negative for
% all nonzero alpha, but we're going backwards on it so the resulting
% locomotion is the same: a positive x body velocity. I think this result
% is much more apparent from the V-style slope plot than the X-style slope
% plot. For example, you can see at a glance that a gait that always
% remains in the negative alpha region will have a negative-theta,
% positive-x body velocity.
%%%%%%%%%%%%%%%%%%%

A_sparse = get_A_of_alpha(s, a_grid(1:2:31));
slope_plot(a_grid(1:2:31),A_sparse,false)
slope_plot(a_grid(1:2:31),A_sparse,true)

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
    sol = asym_solve_gait(s, gait);
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
    sol = asym_solve_gait(s, gait);
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
sol = asym_solve_gait(s, gait);
trajectory_plot(sol,3);
end

hold off
axis equal

%%

% comparing trajcetories of same gait, different phase
figure()
hold on

idx_a = 9; idx_b = 23; % smaller gait
%idx_a = 4; idx_b = 28; % bigger

point_a = a_grid(idx_a)
point_b = a_grid(idx_b)

for phase = 0:pi/16:2*pi
%phase=0;
    
    gait = generate_1D_gait(point_b, 0, phase);
    sol = asym_solve_gait(s, gait);
    traj = trajectory_plot(sol,3);
    
    phase;
    displacement = traj(:,size(traj,2));
    
end

Ax_dif = A.difference(:,1);
Ay_dif = A.difference(:,2);
At_dif = A.difference(:,3);

da = a_grid(2) - a_grid(1);

dx = trapz(Ax_dif(idx_a:idx_b))*da;
dy = trapz(Ay_dif(idx_a:idx_b))*da;
dth = trapz(At_dif(idx_a:idx_b))*da;

for i = 1:3
    g_cerc_vector = [dx dy dth]*i;
    g_circ = [0 -g_cerc_vector(3) g_cerc_vector(1); ...
              g_cerc_vector(3) 0 g_cerc_vector(2); ...
              0 0 0];
    g_exp = expm(g_circ);
    Adif_prediction = mat_to_vec_SE2(g_exp);

    scatter([Adif_prediction(1)],[Adif_prediction(2)],50,'filled')
end

hold off
axis equal

%% plot each of the gaits in time
figure()
hold on

idx_a = 9; idx_b = 23; % smaller gait
%idx_a = 4; idx_b = 28; % bigger

point_a = a_grid(idx_a)
point_b = a_grid(idx_b)

for phase = 0:pi/16:2*pi
    gait = generate_1D_gait(point_b, 0, phase);
    if phase == 2*pi
        plot(linspace(0,3*2*pi),gait{1}(linspace(0,3*2*pi)),'--k')
    elseif phase == pi/2
        plot(linspace(0,3*2*pi),gait{1}(linspace(0,3*2*pi)),'r')
    else
        plot(linspace(0,3*2*pi),gait{1}(linspace(0,3*2*pi)),'c')
    end
end

%idx_a = 9; idx_b = 23; % smaller gait
idx_a = 4; idx_b = 28; % bigger

point_a = a_grid(idx_a)
point_b = a_grid(idx_b)

for phase = 0:pi/16:2*pi
    gait = generate_1D_gait(point_b, 0, phase);
    if phase == 2*pi
        plot(linspace(0,3*2*pi),gait{1}(linspace(0,3*2*pi)),'--k')
    elseif phase == pi/2
        plot(linspace(0,3*2*pi),gait{1}(linspace(0,3*2*pi)),'r')
    else
        plot(linspace(0,3*2*pi),gait{1}(linspace(0,3*2*pi)),'c')
    end
end

set(gca,'XTick',0:pi:6*pi)
hold off

%% Experimenting with integrating Adif to get a prediction of the net
% displacement produced by a gait.
% TODO: explain this process, and maybe a little about the observed error
% between Adif and ODE methods.
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
sol = asym_solve_gait(s, gait);

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
sol = asym_solve_gait(s, gait);

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

Adif_gait_prediction(s, 1, 0)
% this function is unfinished