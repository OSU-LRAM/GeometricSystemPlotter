% creates a band plot of BVI, cBVI, third order contributions to
% displacement

% sanitize
assert(exist('p', 'var'));
assert(exist('s', 'var'));

amplitudes = .25:.5:2.25;
coordinate = '';
is_square = false;

% calculate third-order
[p.to, p.to_opt] = calc_tlb_thirdorder(s,p,is_square);

% generate displacements
truth = cell2mat(cellfun(@(shch) shch.(['G' coordinate])(end,:)', p.G_locus_full, 'UniformOutput', false));
bvi = cell2mat(cellfun(@(shch) shch.(['bvi' coordinate])(end,:)', p.G_locus_full, 'UniformOutput', false));
cbvi = cell2mat(p.(['cBVI' coordinate]));
cbvi_to = cbvi + cell2mat(cellfun(@(shch) shch{1}, p.(['to' coordinate]), 'UniformOutput', false)');

% plotting quantities
truth_points = ones(size(amplitudes));
bvi_points = abs(vecnorm(truth) - vecnorm(truth-bvi))./vecnorm(truth);
cbvi_points = abs(vecnorm(truth) - vecnorm(truth-cbvi))./vecnorm(truth);
cbvi_to_points = abs(vecnorm(truth) - vecnorm(truth-cbvi_to))./vecnorm(truth);

% plot contributions
figure(2);
clf;
hold on;
plot(amplitudes, truth_points, 'k', 'LineWidth', 2);
patch([amplitudes, fliplr(amplitudes)],...
      [bvi_points,... %upper
       zeros(size(amplitudes))],... %lower
      'r', 'FaceAlpha', 0.3);
patch([amplitudes, fliplr(amplitudes)],...
      [cbvi_points,... %upper
       fliplr(bvi_points)],... %lower
      'g', 'FaceAlpha', 0.3);
patch([amplitudes, fliplr(amplitudes)],...
      [cbvi_to_points,... %upper
       fliplr(cbvi_points)],... %lower
      'b', 'FaceAlpha', 0.1);
patch([amplitudes, fliplr(amplitudes)],...
      [truth_points,... %upper
       fliplr(cbvi_to_points)],... %lower
      'b', 'FaceAlpha', 0);
legend('Ground Truth', 'BVI', 'cBVI', 'cBVI with Third Order', 'Error');
title('Contribution of BCH Terms to Displacement Approximation');
xlabel('Amplitude');
ylabel('Contribution %');
[num, dem] = rat(amplitudes);
labels = arrayfun(@(n,d) [num2str(n) '/' num2str(d)], num, dem, 'UniformOutput', false);
xticks(amplitudes);
xticklabels(labels);
yticklabels((100*yticks));
hold off;