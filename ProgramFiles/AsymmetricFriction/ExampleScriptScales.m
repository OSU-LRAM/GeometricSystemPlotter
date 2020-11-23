%% definitions
%%%%%%%%%%%%%%%%%%%

direction_names = ["FF","BF","FB","BB"];

s = sysf_two_link_lowRe;
nlinks = length(s.geometry.linklengths);

s = ensure_connection_and_metric(s); % required for create_grids
s = create_grids(s)
a_grid = s.grid.eval{1};
adot_grid = s.grid.eval{1}*0.1;

%% figure out which non-scaled systems are consistent in which regions of
% the shape-shapechange space
%%%%%%%%%%%%%%%%%%%

% find which system is consistent when
[system_map, count] = identify_piecewise_system(s, adot_grid);

% look at symmetry and missing/redundant areas
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

