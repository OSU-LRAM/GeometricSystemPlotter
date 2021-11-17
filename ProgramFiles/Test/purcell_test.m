%% optimize_so3 test using sysplotter and purcell swimmer
% requirements:
    % sysplotter 'ProgramFiles' directory on matlab path
    % config 'sysplotter_config.mat' file at same level as run directory

%% load sysf
system = 'sysf_three_link_lowRe';
%system = 'sysf_twin_rotor';
sys_calcsystem('calculate', system);
load([system '_calc.mat'], 's');

%% pull out grid, original LC, optimal LC (use 'eval' coordinates)
grid_true = s.grid.eval;
grid_points = cellfun(@(c) unique(c), grid_true, 'UniformOutput', false);
A_orig_full = s.vecfield.eval.content.Avec;
A_opt_full = s.vecfield.eval.content.Avec_optimized;

%% get \theta LC fields for both orig, opt (use as X rotation dim)
% augment with zero fields for Y, Z rotation dim
z = {zeros(length(grid_points{1}))}; %zero field as cell
A_orig = [A_orig_full(3,:); z z; z z];
A_opt_true = [A_opt_full(3,:); z z; z z];

%% use optimize_so3.m on A_orig, grid to produce A_opt_test
% optimize
[grid_test, X, Y, Z, A_opt_test] = optimize_so3(grid_points, A_orig);
% get CCF
DA_orig = calc_ccf(grid_points, A_orig, @rot_mat, @rot_vec);
DA_opt = calc_ccf(grid_points, A_opt_test, @rot_mat, @rot_vec);

%% plot X fields for A_opt (use sysplotter fn for coloring?)
figure(1); clf;
subplot(1,3,1);
quiver(grid_true{1}, grid_true{2}, A_orig{1,1}, A_orig{1,2}, 'Color', [0 0 0]);
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_\theta, original coordinates');
subplot(1,3,2);
quiver(grid_true{1}, grid_true{2}, A_opt_true{1,1}, A_opt_true{1,2}, 'Color', [0 0 0]);
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_\theta, true minimum perturbation coordinates');
subplot(1,3,3);
quiver(grid_test{1}, grid_test{2}, A_opt_test{1,1}, A_opt_test{1,2}, 'Color', [0 0 0]);
axis('square');
xlim([-3 3]);
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_\theta, test minimum perturbation coordinates');

%% plot CCF in orig, opt
figure(2);clf;
subplot(1,2,1);
contour(grid_test{1}, grid_test{2}, DA_orig{1});
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('Original CCF');
subplot(1,2,2);
contour(grid_test{1}, grid_test{2}, DA_opt{1});
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('Optimal CCF');