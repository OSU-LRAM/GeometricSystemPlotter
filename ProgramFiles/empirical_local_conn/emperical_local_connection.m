% Empirical Local Connection calculation
%
% standalone meant to be used after a data cleaner. 
% 
% x_pos, y_pos ----------- [(n timesteps) x (m links)] array of link positions
% time_array ------------- [t] nx1 array of timestamp reported by mocap
% grid_axis -------------- [xmin xmax ymin ymax] boundaries for shapespace evaluation
% num_grid_points -------- scalar. will create pxp square grid of points
% downsample_rate -------- integer >= 1. for selecting only every few points.
%                          Arrays will be re-defined using (1:downsample:end)
% debug ------------------ [plot_neighborhood,plot_valid_connection] 
%                          Defaults to 0. 
%                          Plots the coordinates in shape space of points 
%                          within neighborhood of current grid point, and 
%                          separately if the connection value is valid. 
%                          only for debugging! Very laggy. 

function output = emperical_local_connection(x_pos,y_pos,time_array,grid_axis,num_grid_points,save_file_name,downsample_rate,debug)

% extract alphas, alpha_dots, and body velocity from mocap data
[a1,a2,a1d,a2d,xbd,ybd,tbd] = threeLink_info_from_mocap(x_pos,y_pos,time_array,downsample_rate);

% define the range_radius based on the grid size. This can change to any value.
range_radius = min([(grid_axis(2)-grid_axis(1))*(0.9/(num_grid_points-1)),(grid_axis(4)-grid_axis(3))*(0.9/(num_grid_points-1))]);

% Calculate local connection matrices in a form sysplotter can use
[Ax1,Ax2,Ay1,Ay2,Atheta1,Atheta2,pt1,pt2] = A_from_data(a1,a2,a1d,a2d,xbd,ybd,tbd,time_array,grid_axis,num_grid_pts,range_radius,debug(1),debug(2)); %#ok<*ASGLU>

% use point arrays as alphas for sysplotter. Flip for meshgrid instead of ndgrid. 
alpha1 = pt1'; %#ok<*NASGU>
alpha2 = pt2';

% save the local connection in the systems directory as a mat file
save(fullfile([sysplotter_inputpath '\Systems'],[save_file_name,'.mat']),'alpha1','alpha2','Atheta1','Atheta2','Ax1','Ax2','Ay1','Ay2')

% output the file name so it's easier to find
output = [save_file_name,'.mat'];

end