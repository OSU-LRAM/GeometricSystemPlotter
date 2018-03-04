% get alpha and body velocity from mocap data assuming 3-link geometry

function [a1,a2,a1d,a2d,xbd,ybd,tbd] = threeLink_info_from_mocap(x_pos,y_pos,time_array,downsample_rate)

% set up all the zeros we'll need
a1 = NaN*ones(size(t_span,1),1);
a2 = NaN*ones(size(t_span,1),1);

xb = NaN*ones(size(t_span,1),1);
yb = NaN*ones(size(t_span,1),1);
tb = NaN*ones(size(t_span,1),1);

% use existing position data to extract body and shape data
for i = 1:size(x_pos,1)
    
    % body veloticy assumes four markers and three links   
    xb(i) = (x_pos(i,2)+x_pos(i,3))/2;
    yb(i) = (y_pos(i,2)+y_pos(i,3))/2;
    tb(i) = atan2(y_pos(i,3)-y_pos(i,2),x_pos(i,3)-x_pos(i,2));  
    
    % shape velocity adjusts for body angle
    a1(i) = atan2((y_pos(i,4)-y_pos(i,3)),(x_pos(i,4)-x_pos(i,3))) - tb(i);
    a2(i) = atan2((y_pos(i,1)-y_pos(i,2)),(x_pos(i,1)-x_pos(i,2))) - tb(i);
    
end

% fill in any missing mocap data

xb = fillmissing(xb,'linear');
yb = fillmissing(yb,'linear');
tb = fillmissing(tb,'linear');
a1 = fillmissing(a1,'linear');
a2 = fillmissing(a2,'linear');

% differentiate body and shape arrays to get velocities

xbd = gradient(xb,time_array(2)-time_array(1));
ybd = gradient(yb,time_array(2)-time_array(1));
tbd = gradient(tb,time_array(2)-time_array(1));

a1d = gradient(a1,time_array(2)-time_array(1));
a2d = gradient(a2,time_array(2)-time_array(1));

% downsample the data to make it smoother

% xb = xb(1:downsample_rate:end); % include if you want the body position data out
% yb = yb(1:downsample_rate:end);
% tb = tb(1:downsample_rate:end);

a1 = a1(1:downsample_rate:end);
a2 = a2(1:downsample_rate:end);

xbd = xbd(1:downsample_rate:end);
ybd = ybd(1:downsample_rate:end);
tbd = tbd(1:downsample_rate:end);

a1d = a1d(1:downsample_rate:end);
a2d = a2d(1:downsample_rate:end);

end
