function [alpha,t] = waypoints_along_gait(gaitname,n)
% Get a set of n waypoints along the gait. Gaitname should be the filename
% for the gait, including the shchf_ heading, and the .mat extension can be
% included or neglected without affecting the function of the program

% For now, this program assumes that there is only a single gait, with one
% segment (and works with segment 1 of gait 1 if more are defined in the
% file).

    % Get the paths for the sysplotter directories
    load sysplotter_config
    
    % Load the .mat file for the gait
    load(fullfile(datapath,gaitname))
    
    % Generate n points evenly spaced along the gait
    time_range = p.time_def{1}{1};
    t = linspace(time_range(1), time_range(2),n);
    
    % Find the gait points at these times
    alpha = p.phi_def{1}{1}(t);
    
end