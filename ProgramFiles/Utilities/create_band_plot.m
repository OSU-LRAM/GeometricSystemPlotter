% plots cBVI w.r.t. amplitude, with a third-order bounding band

%% data, styling (change me)
A=0.25:0.25:2;
phi=pi/4:pi/4:2*pi;
coords = '_opt';

cbvi.Color = [0 0 0];
cbvi.LineWidth = 2;
cbvi.LineStyle = '-';

to.FaceColor = [1 1 1]*0.4;
to.FaceAlpha = 0.8;
to.LineStyle = 'none';

bound.FaceColor = [234 14 30]/255;
bound.FaceAlpha = 0.5;
bound.LineStyle = 'none';

%% data generation, plotting

% sanitize
assert(exist('s', 'var'));
assert(exist('p', 'var'));

% get clean cbvi data
cbvi_data = cellfun(@norm, p.(['cBVI' coords]));
cbvi_data = reshape(cbvi_data, [length(phi), length(A)]);
cbvi_plot_data = cbvi_data(1,:);

% compute third order effects
[p.to, p.to_opt] = calc_tlb_thirdorder(s,p,false);
to_data = cellfun(@(to_struct) norm(to_struct{1}), p.(['to' coords]));
to_data = reshape(to_data, [length(phi), length(A)]);
to_true_data = max(to_data, [], 1);
poly_x = [A fliplr(A)];
poly_y = [cbvi_plot_data + to_true_data, fliplr(cbvi_plot_data - to_true_data)];
to_plot_data = polyshape(poly_x, poly_y);

% generate bound
[~, cBVI_fun, to_fun] = bound_third_order(s, [0 0], @A_est_center, @cBVI_est_taylor, 0.2);
to_est_data = vecnorm(cell2mat(arrayfun(to_fun, 2*A, 'UniformOutput', false)));
cbvi_est_data = vecnorm(cBVI_fun(2*A));
to_est_corr = to_est_data./cbvi_est_data.*cbvi_plot_data;
poly_y = [cbvi_est_data + to_est_data, fliplr(cbvi_est_data - to_est_data)];
bound_data = polyshape(poly_x, poly_y);

% plotting
figure(1);
clf;
hold on;
hc = plot(A, cbvi_plot_data, 'Color', cbvi.Color, 'LineWidth', cbvi.LineWidth, 'LineStyle', cbvi.LineStyle);
he = plot(A, cbvi_est_data, 'Color', bound.FaceColor, 'LineWidth', cbvi.LineWidth, 'LineStyle', cbvi.LineStyle);
ht = plot(to_plot_data, 'FaceColor', to.FaceColor, 'LineStyle', to.LineStyle, 'FaceAlpha', to.FaceAlpha);
hb = plot(bound_data, 'FaceColor', bound.FaceColor, 'LineStyle', bound.LineStyle, 'FaceAlpha', bound.FaceAlpha);
set(gca, 'Children', flipud(get(gca, 'Children'))); %flip drawing order
hold off;
xlabel('Gait Amplitude');
ylabel('Resulting Displacement Norm');
legend([hc ht he hb], 'cBVI Estimate', 'Maximum Third Order Contribution', 'cBVI Taylor Series Approximation', 'Estimated Third Order Bound');
[num, dem] = rat(A);
labels = arrayfun(@(n,d) [num2str(n) '/' num2str(d)], num, dem, 'UniformOutput', false);
xticks(A);
xticklabels(labels);
xlim([min(A) max(A)]);
ylim([0 1.2*max(poly_y)]);