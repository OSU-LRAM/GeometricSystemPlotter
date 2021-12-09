direction_names = ["FF","BF","FB","BB"];

s = sysf_two_link_lowRe;
nlinks = length(s.geometry.linklengths);

% FOR NULL CASE
% s.physics.drag_bw_ratio = 1;

s = ensure_connection_and_metric(s); % required for create_grids
s = create_grids(s)
a_grid = s.grid.eval{1};
adot_grid = s.grid.eval{1}*0.1;

    
%% Demonstrate getting velocities and plotting as vectors in body and link
% origins
%%%%%%%%%%%%%%%%%%%
a = 1.8;
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
velocity_diagram(s, lvel, bvel, a, adot, 0, 1, 1);
axis off

a = -1.8;
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
velocity_diagram(s, lvel, bvel, a, adot, 0, 1, 1);
axis off

a = 1.8;
adot = -0.5;
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
velocity_diagram(s, lvel, bvel, a, adot, 0, 1, 1);
axis off

a = -1.8;
adot = -0.5;
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
velocity_diagram(s, lvel, bvel, a, adot, 0, 1, 1);
axis off

%% Axis Labels
display.aspect_ratio = 0.05;
display.sharpness = 0.1;

figure()
axis manual
ax = gca;
ax.NextPlot = 'replaceChildren';
xlim([-1, 1]);
ylim([-1, 1]);

a = 0.5;
B = fat_chain(s.geometry, a, display);
patch(B(1,:),B(2,:),[1 0 0])
axis off

%% Some example gaits.

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

% Get A as a function of alpha (split into A.positive and A.negative
A = get_A_of_alpha(s, a_grid);

%% Make a plot for each of x, y, and theta components of A. Plot A+ and A-

% directional glyphs
figure()
title("A- and A+");
tiledlayout(3,1)

% X plot
ax1 = nexttile;
plot(ax1,a_grid, A.positive(:,1), 'r>', a_grid, A.negative(:,1), 'k<');
ylabel(ax1,'A_x')
xlabel(ax1,'\alpha')
xticks([-2.5:1:2.5])

% Y plot
ax1 = nexttile;
plot(ax1,a_grid, A.positive(:,2), 'r>', a_grid, A.negative(:,2), 'k<');
ylabel(ax1,'A_y')
xlabel(ax1,'\alpha')
xticks([-2.5:1:2.5])

% Theta plot
ax2 = nexttile;
plot(ax2,a_grid, A.positive(:,3), 'r>', a_grid, A.negative(:,3), 'k<');
ylabel(ax2,'A_\theta')
xlabel(ax2, '\alpha')
xticks([-2.5:1:2.5])

figure()
title("Adif");
tiledlayout(2,1)

% X plot
ax1 = nexttile;
plot(ax1,a_grid, A.difference(:,1),'k')
ylabel(ax1,'A_x')
xlabel(ax1, '\alpha')
xticks([-2.5:1:2.5])

% Theta plot
ax2 = nexttile;
plot(ax2,a_grid, A.difference(:,3),'k')
ylabel(ax2,'A_\theta')
xlabel(ax2, '\alpha')
xticks([-2.5:1:2.5])

%% Make slope plots in both styles.
%%%%%%%%%%%%%%%%%%%

A_sparse = get_A_of_alpha(s, a_grid(1:2:31));
% slope_plot(a_grid(1:2:31),A_sparse,false)
% xticks([-2.5:1:2.5])
slope_plot(a_grid(1:2:31),A_sparse,true)

%% Gaits for paper figure section 4

centered_gait = generate_1D_gait(1, 0, 0);
offset_gait = generate_1D_gait(1, 1, 0);
gait = offset_gait;

% Centered Gait
point_b = gait{1}(pi/2);
point_a = gait{1}(3*pi/2);
da = (point_b-point_a)/16;

% New Adiff going from point a to point b
Adiff = get_A_of_alpha(s,point_a:da:point_b).difference;

Ax_dif = Adiff(:,1);
Ay_dif = Adiff(:,2);
At_dif = Adiff(:,3);

dx = trapz(Ax_dif)*da;
dy = trapz(Ay_dif)*da;
dth = trapz(At_dif)*da;

g_circ_vector = [dx dy dth];
g_circ = [0 -g_circ_vector(3) g_circ_vector(1); ...
          g_circ_vector(3) 0 g_circ_vector(2); ...
          0 0 0];
Adif_prediction = mat_to_vec_SE2(expm(g_circ))

% compare to a gait that goes from point_a to point_b:
sol = asym_solve_gait(s, gait);

displacement = deval(sol, 2*pi)'

error = Adif_prediction - displacement;
proport_error = error ./ displacement;

figure()
hold on
trajectory_plot(sol,1)
quiver(0,0,1,0,0.01)
x=Adif_prediction(1); y=Adif_prediction(2); th=Adif_prediction(3);
quiver(x,y,x+cos(th),y+sin(th),0.01)
hold off
axis equal

%% Experimenting with integrating Adif to get a prediction of the net
% displacement produced by a gait.
% TODO: explain this process, and maybe a little about the observed error
% between Adif and ODE methods.
%%%%%%%%%%%%%%%%%%%

%% Using the function that packages this calculation into something less
% cumbersome
%%%%%%%%%%%%%%%%%%%

%Adif_gait_prediction(s, 1, 0, calc_density)
% this function is unfinished