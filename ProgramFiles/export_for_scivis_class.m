function export_for_scivis_class(core_path, system_name)
%EXPORT_FOR_SCIVIS_CLASS Summary of this function goes here
% Run sysplotter first (and then close it) in order to get the
% sysplotter_inputpath variable if you want it.

% Systems must already have been calculated and cached in Sysplotter for
% this to run.

file = strcat(core_path, '/sysplotter_data/', system_name, '_calc.mat');
load(file, 's');

output_base_name = strcat(core_path, '/sysplotter_data/', system_name, '/', system_name, '_');

fldnm = 'DA_optimized';
idx_names = ['x','y','theta'];
for idx = 1:3
    csvwrite(strcat(output_base_name, fldnm, '_', idx_names(idx), '.csv'), s.(fldnm){idx});
end

fldnm = 'DA';
for idx = 1:3
    csvwrite(strcat(output_base_name, fldnm, '_', idx_names(idx), '.csv'), s.(fldnm){idx});
end

idx_names = ['a', 'b'];
for idx = 1:2
    csvwrite(strcat(output_base_name, 'grid', '_', idx_names(idx), '.csv'), s.grid.eval{idx});
end

f = figure('visible','off');
contour(s.grid.eval{1},s.grid.eval{2},s.DA_optimized{1});
axis([s.grid_range]);
axis equal;
xticks([-2.5 -1.5 -0.5 0.5 1.5 2.5]);
yticks([-2.5 -1.5 -0.5 0.5 1.5 2.5]);
print(strcat(output_base_name,'_sample_plot.png'),'-dpng');

end