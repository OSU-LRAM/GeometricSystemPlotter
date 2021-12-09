function [consistent_with_SSS, consistent_subsystem, trajectory, t_solution] = when_do_scales_change_things(s,gait, gait_repeats, frames_per_cycle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nlinks = length(s.geometry.linklengths);

sol = asym_solve_gait(s, gait);

% when is it going forward all links

directions_possible = ff2n(nlinks);
i = 1; % SSS
links_backwards = directions_possible(i,:);
friction_directions = ones(1,nlinks) - 2 * links_backwards;

gait_period = 2*pi;
nframes = gait_repeats * frames_per_cycle;

t_solution = linspace(0, gait_period*gait_repeats, nframes);
g_solution = deval(sol, mod(t_solution, gait_period));

gait_displacement = vec_to_mat_SE2(deval(sol, gait_period));

trajectory = zeros(size(g_solution));
consistent_with_SSS = zeros(1,nframes);
consistent_subsystem = zeros(1,nframes);

for frame = 1:nframes
    t = t_solution(frame);
    gaits_finished = floor(t/gait_period);
    
    accumulated_displacement = gait_displacement ^ gaits_finished;
    g = accumulated_displacement * vec_to_mat_SE2(g_solution(:, frame));
    
    
    trajectory(:, frame) = mat_to_vec_SE2(g);
    
    % now color code with the current consistency
    shape = gait{1}(t);
    [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
    shapechange = gait{2}(t);
    link_directions = [0 0 0];
    [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry,s.physics, shape);
    body_vel = -A * shapechange;
    for link = 1:length(s.geometry.linklengths)
        link_body_vel = J_full{link} * [body_vel; shapechange];
        dx_direction = sign(link_body_vel(1));
        link_directions(link) = dx_direction;
    end
    all_links_agree = all(~link_directions | ~(link_directions - friction_directions));
    if all(link_directions)
        % get the numbered subsystem currently being observed according to
        % link motion
        consistent_subsystem(frame) = bin2dec(num2str( (-link_directions+1)/2 )) + 1;
    else
        % if any link direction is currently zero, use a dummy value
        consistent_subsystem(frame) = 9;
    end
    
    consistent_with_SSS(frame) = all_links_agree;
    
end

end

