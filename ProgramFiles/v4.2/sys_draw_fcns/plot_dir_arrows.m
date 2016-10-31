function arrow_lines = plot_dir_arrows(x,y,n,varargin)
%Plots the direction arrows for a parametric path [x(t),y(t)] at n points

if n ~=0

    %ensure input vectors are columns
    if (size(x,1)<=size(x,2))
        x = x';
        y = y';
    end

    %figure out size of curve, to scale arrowheads
    x_span = max(x) - min(x);
    y_span = max(y) - min(y);
    max_span = max([x_span y_span]);


    %use number of arrows to draw to determine indices of arrow locations
    approx_locs = length(x)*(1/(2*n)+(0:n-1)/n);

    %intergerize the approximate locations, aiming to go slightly ahead of the
    %ideal index for the tip
    tip_locs = ceil(approx_locs)+1;
    %index one back from tip location, to get tail of tangent vector
    tail_locs = tip_locs-1;

    %convert from indices to values for tip and tail positions
    tip_pos = [x(tip_locs) y(tip_locs)];
    tail_pos = [x(tail_locs) y(tail_locs)];

    %get vector from tip to tail (negative tangent vector)
    tan_vec = tail_pos-tip_pos;

    %get magnitude of each negative tangent vector
    mag_tan_vec = hypot(tan_vec(:,1),tan_vec(:,2));

    %scale the negative tangent vector to 1/8 of the maximum span of the curve
    tan_vec_scaled = tan_vec  ./ (mag_tan_vec * [1 1]) * max_span/8;

    % find the outboard spacing of the arrowheads, using a 1.6 golden rule
    % ratio
    normal_comp = [-tan_vec_scaled(:,2) tan_vec_scaled(:,1)]*1/1.6;

    %find the postition of the two barb ends
    first_barb_end = tip_pos+tan_vec_scaled+normal_comp;
    second_barb_end = tip_pos+tan_vec_scaled-normal_comp;

    %create 2-point lines from tip to barb ends
    barbs_x = [ tip_pos(:,1) first_barb_end(:,1);
                tip_pos(:,1) second_barb_end(:,1)];

    barbs_y = [ tip_pos(:,2) first_barb_end(:,2);
                tip_pos(:,2) second_barb_end(:,2)];

    %plot the arrowheads
    arrow_lines = line(barbs_x',barbs_y',varargin{:});
    
end
