%% styling (change me!)
is_opt = false;
is_square = true;

net.LineStyle = 'none';
net.LineWidth = 1;
net.Color = [0 0 0];
net.Marker = '.';
net.MarkerSize = 15;
net.MarkerFaceColor = [0 0 0];

cbvi.LineStyle = 'none';
cbvi.LineWidth = 1;
cbvi.Color = [0 0 0];
cbvi.Marker = 'o';
cbvi.MarkerSize = 10;

to.LineStyle = 'none';
to.LineWidth = 2;
to.Color = [234 14 30]/255;
to.Marker = 'x';
to.MarkerSize = 10;

%% sanitize
assert(exist('s','var'));
assert(exist('p','var'));

%% plot restyling
% get open axes
axes = gca;

% add third-order effects
add_tlb_plot(s,p,is_opt,is_square,true);

% store net, third order
to_data = [];

% use DisplayName to sort data
for i = 1:length(axes.Children)
    child = axes.Children(i);
    if contains(child.DisplayName, 'Third Order')
        % third-order effects
        edit_properties(child, to);
        to_data = [to_data [child.XData; child.YData]];
    elseif contains(child.DisplayName, 'cBVI')
        % cbvi
        edit_properties(child, cbvi);
    elseif contains(child.DisplayName, 'Net')
        % ground truth
        edit_properties(child, net);
        net_data = [child.XData; child.YData];
    end
end

% add corresp. lines
empties = find(cellfun(@isempty, p.cBVI));
net_data(:,empties) = []; %clean displacements without cBVI's
to_data = fliplr(to_data); %flipped for some reason
for i = 1:length(net_data)
    x = [net_data(1,i) to_data(1,i)];
    y = [net_data(2,i) to_data(2,i)];
    line(axes, x, y, 'Color', 'k', 'LineStyle', ':');
end

% add bound
[len, cBVI_fun, to_fun] = bound_third_order(s, [0 0],...
                                            @(a,b) A_est_center(a,b,is_opt),...
                                            @(a,b) cBVI_est_taylor(a,b,is_opt),...
                                            0.5, is_square);
amps = linspace(0, len, 20);
bound_data = -cell2mat(arrayfun(@(a) cBVI_fun(a), amps, 'UniformOutput', false)) +...
             cell2mat(arrayfun(@(a) to_fun(a), amps, 'UniformOutput', false));
bound_data_neg = bound_data.*[1 -1 1]';
%line(axes, bound_data(1,:), bound_data(2,:), 'Color', [234 14 30]/255, 'LineStyle', '--');
%line(axes, bound_data_neg(1,:), bound_data_neg(2,:), 'Color', [234 14 30]/255, 'LineStyle', '--');

function edit_properties(child, type)
    fields = fieldnames(type);
    for i = 1:length(fields)
        field = fields{i};
        child.(field) = type.(field);
    end
end