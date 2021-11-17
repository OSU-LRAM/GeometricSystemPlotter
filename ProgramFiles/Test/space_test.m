%% optimize_so3 test on a reaction wheel satellite
% create samples in ea. shape dimension
num_samples=31; bound=pi/4;
dim = linspace(-bound,bound,num_samples);
samples = {dim,dim};
% create local connection
z = zeros(num_samples);
o = ones(num_samples);
% the following are valid positive configurations
A_orig = {z, z; -o, z; z, -o};
%A_orig = {-o z; z z; z -o};
%A_orig = {-o z; z -o; z z};
% optimize
[grid, X, Y, Z, A_opt] = optimize_so3(samples, A_orig);

%% compare lengths
X_orig = [A_orig{1,1}(:)';A_orig{1,2}(:)'];
Y_orig = [A_orig{2,1}(:)';A_orig{2,2}(:)'];
Z_orig = [A_orig{3,1}(:)';A_orig{3,2}(:)'];
n_orig = norm([mean(vecnorm(X_orig)) mean(vecnorm(Y_orig)) mean(vecnorm(Z_orig))]);
disp(['Original norm-average metric: ', num2str(n_orig)]);

X_opt = [A_opt{1,1}(:)';A_opt{1,2}(:)'];
Y_opt = [A_opt{2,1}(:)';A_opt{2,2}(:)'];
Z_opt = [A_opt{3,1}(:)';A_opt{3,2}(:)'];
n_opt = norm([mean(vecnorm(X_opt)) mean(vecnorm(Y_opt)) mean(vecnorm(Z_opt))]);
disp(['Optimal norm-average metric: ', num2str(n_opt)]);

%% test gaits
% create gait
t_vec = linspace(0, 2*pi);
alpha = @(t) 0.5*bound*[cos(t); sin(t)];
d_alpha = @(t) 0.5*bound*[-sin(t); cos(t)];
% vector for plotting
alpha_vec = alpha(t_vec);
% integrate over connection
[g_circ_orig, g_orig] = so3_integrator(t_vec, alpha, d_alpha, grid, A_orig);
[g_circ_opt, g_opt] = so3_integrator(t_vec, alpha, d_alpha, grid, A_opt);

%% constraint curvature
DA_orig = calc_ccf(samples, A_orig, @rot_mat, @rot_vec);
DA_opt = calc_ccf(samples, A_opt, @rot_mat, @rot_vec);

%% plot connections
figure(1);clf;

subplot(2,3,1);
hold on;
quiver(grid{1},grid{2},A_orig{1,1},A_orig{1,2}, 'Color', [0 0 0]);
plot(alpha_vec(1,:), alpha_vec(2,:), 'Color', [234 14 30]/255);
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_x, original coordinates');
subplot(2,3,2);
hold on;
quiver(grid{1},grid{2},A_orig{2,1},A_orig{2,2}, 'Color', [0 0 0]);
plot(alpha_vec(1,:), alpha_vec(2,:), 'Color', [234 14 30]/255);
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_y, original coordinates');
subplot(2,3,3);
hold on;
quiver(grid{1},grid{2},A_orig{3,1},A_orig{3,2}, 'Color', [0 0 0]);
plot(alpha_vec(1,:), alpha_vec(2,:), 'Color', [234 14 30]/255);
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_z, original coordinates');

subplot(2,3,4);
hold on;
quiver(grid{1},grid{2},A_opt{1,1},A_opt{1,2}, 'Color', [0 0 0]);
plot(alpha_vec(1,:), alpha_vec(2,:), 'Color', [234 14 30]/255);
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_x, optimal coordinates');
subplot(2,3,5);
hold on;
quiver(grid{1},grid{2},A_opt{2,1},A_opt{2,2}, 'Color', [0 0 0]);
plot(alpha_vec(1,:), alpha_vec(2,:), 'Color', [234 14 30]/255);
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_y, optimal coordinates');
subplot(2,3,6);
hold on;
quiver(grid{1},grid{2},A_opt{3,1},A_opt{3,2}, 'Color', [0 0 0]);
plot(alpha_vec(1,:), alpha_vec(2,:), 'Color', [234 14 30]/255);
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_z, optimal coordinates');

%% plot log coords
figure(2); clf;

subplot(2,3,1);
plot(t_vec, g_circ_orig(1,:), 'Color', [234 14 30]/255);
xlabel('Time');
ylabel('Displacement');
title('X Diplacement, Original Coordinates');
subplot(2,3,2);
plot(t_vec, g_circ_orig(2,:), 'Color', [234 14 30]/255);
xlabel('Time');
ylabel('Displacement');
title('Y Diplacement, Original Coordinates');
subplot(2,3,3);
plot(t_vec, g_circ_orig(3,:), 'Color', [234 14 30]/255);
xlabel('Time');
ylabel('Displacement');
title('Z Diplacement, Original Coordinates');

subplot(2,3,4);
plot(t_vec, g_circ_opt(1,:), 'Color', [234 14 30]/255);
xlabel('Time');
ylabel('Displacement');
title('X Diplacement, Optimal Coordinates');
subplot(2,3,5);
plot(t_vec, g_circ_opt(2,:), 'Color', [234 14 30]/255);
xlabel('Time');
ylabel('Displacement');
title('Y Diplacement, Optimal Coordinates');
subplot(2,3,6);
plot(t_vec, g_circ_opt(3,:), 'Color', [234 14 30]/255);
xlabel('Time');
ylabel('Displacement');
title('Z Diplacement, Optimal Coordinates');

%% plot CCF
figure(3); clf;
subplot(2,3,1);
contour(grid{1}, grid{2}, DA_orig{1});
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('X CCF, Original Coordinates');
subplot(2,3,2);
contour(grid{1}, grid{2}, DA_orig{2});
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('Y CCF, Original Coordinates');
subplot(2,3,3);
contour(grid{1}, grid{2}, DA_orig{3});
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('Z CCF, Original Coordinates');

subplot(2,3,4);
contour(grid{1}, grid{2}, DA_opt{1});
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('X CCF, Optimal Coordinates');
subplot(2,3,5);
contour(grid{1}, grid{2}, DA_opt{2});
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('Y CCF, Optimal Coordinates');
subplot(2,3,6);
contour(grid{1}, grid{2}, DA_opt{3});
axis('equal');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('Z CCF, Optimal Coordinates');
