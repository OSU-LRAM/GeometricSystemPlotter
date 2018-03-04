

function [Ax1,Ax2,Ay1,Ay2,Atheta1,Atheta2,pt1,pt2] = A_from_data(alpha_1,alpha_2,alpha_1_dot,alpha_2_dot,x_dot_body,y_dot_body,th_dot_body,time_array,grid_axis,num_grid_pts,range_radius,plot_neighborhood,plot_values)
                                                %    A_from_data(a1,a2,xb,yb,tb,time_array,grid_axis,num_grid_points);
if ~exist('plot_neighborhood','var')
    plot_neighborhood = 0;
end
if ~exist('plot_values','var')
    plot_values = 0;
end
    
a1 = alpha_1;
a2 = alpha_2;
ad1 = alpha_1_dot;
ad2 = alpha_2_dot;

% start analysis given grid parameters
% grid_array = linspace(grid_min,grid_max,num_grid_pts);
% grid_array = linspace(min([min(a1),min(a2)])+0.1,max([max(a1),max(a2)])-0.1,num_grid_pts);
grid_a1_array = linspace(grid_axis(1),grid_axis(2),num_grid_pts);
grid_a2_array = linspace(grid_axis(3),grid_axis(4),num_grid_pts);
[pt1,pt2] = meshgrid(grid_a1_array,grid_a2_array);

% all the zeroes to make Matlab happy
dAx = zeros(size(pt1,1),size(pt2,2));
dAy = zeros(size(pt1,1),size(pt2,2));
dAt = zeros(size(pt1,1),size(pt2,2));
Ax1 = dAx;
Ax2 = dAx;
Ay1 = dAx;
Ay2 = dAx;
Atheta1 = dAx;
Atheta2 = dAx;

% if ~exist('range_radius','var')
%     % make the circles overlap but not double up
%     range_radius = min([(grid_a1_array(2)-grid_a1_array(1))*0.9,(grid_a2_array(2)-grid_a2_array(1))*0.9]);
% end

% at each grid point, calculate the local connection height value
for pti = 1:num_grid_pts
    for ptj = 1:num_grid_pts
        
        pt = [pt1(pti,ptj),pt2(pti,ptj)];
        
        % get all the neighborhood points around current grid point
        [d1,d2] = neighborhood_points(pt,[a1,a2],grid_axis,range_radius,'square',0,plot_neighborhood);
        %neighborhood_points(home_point,outside_points,range_radius,bound_type,edge_adjust,plot_test)

        
        % run the empirical calculation
        Ax_mat_out = empirical_local_conn_calculation(x_dot_body,[ad1,ad2],d1,d2);
        Ay_mat_out = empirical_local_conn_calculation(y_dot_body,[ad1,ad2],d1,d2);
        At_mat_out = empirical_local_conn_calculation(th_dot_body,[ad1,ad2],d1,d2);
        
        A1ad2x = Ax_mat_out(4); % select dA1/da2 and dA2/da1
        A2ad1x = Ax_mat_out(5);
        dAx(pti,ptj) = A2ad1x - A1ad2x; % surface value for curl
        Ax1(pti,ptj) = Ax_mat_out(1);
        Ax2(pti,ptj) = Ax_mat_out(2);
        
        A1ad2y = Ay_mat_out(4);
        A2ad1y = Ay_mat_out(5);
        dAy(pti,ptj) = A2ad1y - A1ad2y;
        Ay1(pti,ptj) = Ay_mat_out(1);
        Ay2(pti,ptj) = Ay_mat_out(2);
        
        A1ad2t = At_mat_out(4);
        A2ad1t = At_mat_out(5);
        dAt(pti,ptj) = A2ad1t - A1ad2t;
        Atheta1(pti,ptj) = At_mat_out(1);
        Atheta2(pti,ptj) = At_mat_out(2);
        
        if plot_values
%            
            figure(15); 
            hold on
            if sum(Ax_mat_out) == 0 %isnan(Ax1(pti,ptj))
%                 subplot(1,2,2); 
%                 hold on
                plot(pt(1),pt(2),'r+','LineWidth',2)
            else
%                 subplot(1,2,2); 
%                 hold on
                plot(pt(1),pt(2),'b+','LineWidth',2)
            end
            title('Ax1: \color{blue}Valid output \color{red}Contains NaN')
            axis equal
            axis([min(min(pt1)) max(max(pt1)) min(min(pt2)) max(max(pt2))])
            
            
                        
            
%             subplot(3,2,2); hold on
%             if isnan(Ax2(pti,ptj))
%                 plot(pt(1),pt(2),'r+','LineWidth',3)
%             else
%                 plot(pt(1),pt(2),'b+','LineWidth',3)
%             end
%             title('Ax2')
%             axis equal
%             axis([min(min(pt1)) max(max(pt1)) min(min(pt2)) max(max(pt2))])
%             
%             subplot(3,2,3); hold on
%             if isnan(Ay1(pti,ptj))
%                 plot(pt(1),pt(2),'r+','LineWidth',3)
%             else
%                 plot(pt(1),pt(2),'b+','LineWidth',3)
%             end
%             title('Ay1')
%             axis equal
%             axis([min(min(pt1)) max(max(pt1)) min(min(pt2)) max(max(pt2))])
%                         
%             
%             subplot(3,2,4); hold on
%             if isnan(Ay2(pti,ptj))
%                 plot(pt(1),pt(2),'r+','LineWidth',3)
%             else
%                 plot(pt(1),pt(2),'b+','LineWidth',3)
%             end
%             title('Ay2')
%             axis equal
%             axis([min(min(pt1)) max(max(pt1)) min(min(pt2)) max(max(pt2))])
%                         
%             subplot(3,2,5); hold on
%             if isnan(Atheta1(pti,ptj))
%                 plot(pt(1),pt(2),'r+','LineWidth',3)
%             else
%                 plot(pt(1),pt(2),'b+','LineWidth',3)
%             end
%             title('A\theta1')
%             axis equal
%             axis([min(min(pt1)) max(max(pt1)) min(min(pt2)) max(max(pt2))])
%                         
%             
%             subplot(3,2,6); hold on
%             if isnan(Atheta2(pti,ptj))
%                 plot(pt(1),pt(2),'r+','LineWidth',3)
%             else
%                 plot(pt(1),pt(2),'b+','LineWidth',3)
%             end
%             title('A\theta2')
%             axis equal
%             axis([min(min(pt1)) max(max(pt1)) min(min(pt2)) max(max(pt2))])
            
%             pause(0.001)
            
        end
        
    end
end


end