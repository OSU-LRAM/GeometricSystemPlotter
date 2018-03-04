% Empirical Local Connection Data Cleaner
%
% Reads in motion capture data and modifies the axis and point order so the
% empirical local connection calculations work correctly. 
%
% MOCAP DATA MUST BE CONVERTED TO THE FOLLOWING FORMAT.
% x,y,z ---- (n x m) arrays where n is the number of timesteps and m is the
%            number of markers
% t -------- (n x 1) time array


function [new_x,new_y,new_z,t] = clean_mocap_data(x,y,z,t,fix_axis_order,fix_point_order)

num_markers = size(x,2);

% fill in any missing data. 
xpts = fillmissing(x,'linear',1);
ypts = fillmissing(y,'linear',1);
zpts = fillmissing(z,'linear',1);

% xpts = NaN*ones(size(x,1),num_markers);
% ypts = NaN*ones(size(y,1),num_markers);
% zpts = NaN*ones(size(z,1),num_markers);

points_ok = 0;
axis_ok = 0;

% Re-order the axes displayed. 
if ~exist('fix_axis_order','var')
    
%     fix_axis_order = [1,2,3];
    disp('Plotting sample configuration for axis evaluation')
    while axis_ok ~= 1
        fig = figure(120); clf;  hold on %#ok<NASGU>
%         plot3(xpts(1,1),ypts(1,1),zpts(1,1),'ro','MarkerSize',4,'LineWidth',3)
        for j = 1:num_markers
            plot3(xpts(1,j),ypts(1,j),zpts(1,j),'k+','MarkerSize',4,'LineWidth',3)
            xlabel('x'); ylabel('y'); zlabel('z');
%             text(xpts(1,j),ypts(1,j),num2str(j))
            axis equal
        end
        grid on
        axis equal
        
        axis_ok = input('Does this match the orientation expected? 0 or 1 --> ');
        
        
        if axis_ok == 1
            continue
            close(fig) %#ok<UNRCH>
        elseif axis_ok == 2
            error('User stopped analysis. Sample layout was too different.')
        else
            fix_axis_order = input('Enter new axis layout. [dim1 dim2 dim3] --> ');
            
            for i = 1:size(xpts,2)
                
                dummy = [xpts(:,i),ypts(:,i),zpts(:,i)];
                
                dummyx(:,i) = dummy(:,fix_axis_order(1));
                dummyy(:,i) = dummy(:,fix_axis_order(2));
                dummyz(:,i) = dummy(:,fix_axis_order(3));
                
            end
            xpts = dummyx;
            ypts = dummyy;
            zpts = dummyz;
            
        end
    end
else
    for i = 1:size(xpts,2)
        dummy = [xpts(:,i),ypts(:,i),zpts(:,i)];
                
                dummyx(:,i) = dummy(:,fix_axis_order(1));
                dummyy(:,i) = dummy(:,fix_axis_order(2));
                dummyz(:,i) = dummy(:,fix_axis_order(3));
    end
    xpts = dummyx;
    ypts = dummyy;
    zpts = dummyz;    
end

% Re-order the points in space.
if ~exist('fix_point_order','var')
    
%     fix_point_order = 1:num_markers;
    disp('Plotting sample configuration for point order evaluation')
    while points_ok ~= 1
        fig = figure(121); clf;  hold on %#ok<NASGU>
        
        for j = 1:num_markers
            plot(xpts(1,j),ypts(1,j),'k+','MarkerSize',4)
            text(xpts(1,j),ypts(1,j),num2str(j))
            axis equal
        end
        xlabel('x'); ylabel('y'); zlabel('z');
        grid on
        axis equal
        
        points_ok = input('Does this order match the expected? 0 or 1 --> ');
        
        
        if points_ok == 1
            continue
            close(fig) %#ok<UNRCH>
        elseif points_ok == 2
            error('User stopped analysis. Sample layout was too different.')
        else
            fix_point_order = input('Enter order from left to right displayed. [p1 p2 ... pn] --> ');
            
            for i = 1:size(xpts,2)
                
                dummyx(:,i) = xpts(:,fix_point_order(i)); %#ok<*AGROW>
                dummyy(:,i) = ypts(:,fix_point_order(i));
                dummyz(:,i) = zpts(:,fix_point_order(i));
                
            end
            xpts = dummyx;
            ypts = dummyy;
            zpts = dummyz;
            
        end
    end
else
    for i = 1:size(xpts,2)
        dummyx(:,i) = xpts(:,fix_point_order(i)); %#ok<*AGROW>
        dummyy(:,i) = ypts(:,fix_point_order(i));
        dummyz(:,i) = zpts(:,fix_point_order(i));
    end
    xpts = dummyx;
    ypts = dummyy;
    zpts = dummyz;    
end


new_x = xpts;
new_y = ypts;
new_z = zpts;

end