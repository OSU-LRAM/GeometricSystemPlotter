% find the points in the local neighborhood of a given center point.
% output gives two delta matrices. Each has one column defining the index
% in the original outside_points matrix, and a second column defining the
% distance to that point in delta1 or delta2 direction. d1 corresponds to x
% distance and d2 to y.

function [d1,d2] = neighborhood_points(home_point,outside_points,range_radius,bound_type,edge_adjust,plot_test)

n = size(outside_points,1);
a = outside_points;
pt = home_point;
srange = range_radius;

d1 = [];
d2 = [];

if edge_adjust
    
    % set default ranges
    srange1L = srange;
    srange1R = srange;
    srange2D = srange;
    srange2U = srange;
    
    % adjust ranges
    if pt(1)-srange < min(a(:,1))
        srange1L = pt(1)-min(a(:,1));
        srange1R = 1.5*srange-(pt(1)-min(a(:,1)));
    end
    if pt(1)+srange > max(a(:,1))
        srange1L = 1.5*srange-(max(a(:,1))-pt(1));
        srange1R = max(a(:,1))-pt(1);
    end
    
    if pt(2)-srange < min(a(:,2))
        srange2D = pt(2)-min(a(:,2));
        srange2U = 1.5*srange-(pt(2)-min(a(:,2)));
    end
    if pt(2)+srange > max(a(:,2))
        srange2D = 1.5*srange-(max(a(:,2))-pt(2));
        srange2U = max(a(:,2))-pt(2);
    end
else
    srange1L = srange;
    srange1R = srange;
    srange2D = srange;
    srange2U = srange;
end
    

for k = 1:n    
    
    if strcmp(bound_type,'square')
        if pt(1)-a(k,1) <= srange1L && a(k,1)-pt(1) <= srange1R
            if pt(2)-a(k,2) <= srange2D && a(k,2)-pt(2) <= srange2U
                % add [x_index,delta1] to end matrix
                d1 = [d1;k,(a(k,1)-pt(1))]; %#ok<AGROW>
                % add [x_index,delta2] to end matrix
                d2 = [d2;k,(a(k,2)-pt(2))]; %#ok<AGROW>
            end
        end
        
    elseif strcmp(bound_type,'circle')
        if sqrt((pt(1) - a(k,1))^2 + (pt(2) - a(k,2))^2) <= srange
            % add [x_index,delta] to end matrix
            d1 = [d1;k,-(a(k,1)-pt(1))]; %#ok<AGROW>
            d2 = [d2;k,-(a(k,2)-pt(2))]; %#ok<AGROW>
        end
    else
        error('Error: undefined boundary type.')
    end
end


if plot_test

    plotX = [pt(1)-srange1L,pt(1)+srange1R,pt(1)+srange1R,pt(1)-srange1L,pt(1)-srange1L];
    plotY = [pt(2)-srange2D,pt(2)-srange2D,pt(2)+srange2U,pt(2)+srange2U,pt(2)-srange2D];
    
    figure(12); clf; hold on
    plot(pt(1),pt(2),'k+','LineWidth',2)
    plot(a(:,1),a(:,2),'b+')
    plot(plotX,plotY,'r')
    plot(d1(:,2)+pt(1)*ones(size(d1,1)),d2(:,2)+pt(2)*ones(size(d2,1)),'r+')
    axis equal
    pause(0.01)
    
end

end