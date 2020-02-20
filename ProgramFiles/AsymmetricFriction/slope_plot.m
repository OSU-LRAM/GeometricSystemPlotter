function slope_plot(grid,A,style)
%SLOPE_PLOT Summary of this function goes here
%   Detailed explanation goes here

% NaN in plot makes a line break
buffer_ratio = 0.3;

slope_x_pos = zeros([(length(grid)-1)*3, 1]);
slope_x_neg = zeros(size(slope_x_pos));
slope_y_pos = zeros([(length(grid)-1)*3, 3]);
slope_y_neg = zeros(size(slope_y_pos));
for i=1:length(grid)-1
    j = (i - 1) * 3 + 1;
    spacing = grid(i+1) - grid(i);
    buffer = spacing * buffer_ratio;
    da = spacing - buffer;
    if style
    % one way:
        slope_x_pos(j:j+2) = [grid(i), grid(i)+da/2, NaN];
        slope_x_neg(j:j+2) = [grid(i)-da/2, grid(i), NaN];
    else
    % other way:
        slope_x_pos(j:j+2) = [grid(i)-da/2, grid(i)+da/2, NaN];
        slope_x_neg(j:j+2) = [grid(i)-da/2, grid(i)+da/2, NaN];
    end
    for k=1:3
        dx_pos = A.positive(i, k)*da;
        dx_neg = A.negative(i, k)*da;
        if style
        % one way:
            slope_y_pos(j:j+2,k) = [0,  dx_pos/2, NaN];
            slope_y_neg(j:j+2,k) = [-dx_neg/2, 0, NaN];
        else
        % other way:
            slope_y_pos(j:j+2,k) = [-dx_pos/2  dx_pos/2 NaN];
            slope_y_neg(j:j+2,k) = [-dx_neg/2  dx_neg/2 NaN];
        end
    end
end

figure()
tiledlayout(3,1)
for i=1:3
    ax = nexttile;
    plot(ax, slope_x_pos, slope_y_pos(:,i), 'r', slope_x_neg, slope_y_neg(:,i), 'k');
end

end

