function trajectory_plot(solution, gait_cycles)
%TRAJECTORY_PLOT Summary of this function goes here
%   Detailed explanation goes here

gait_period = 2*pi;
frames_per_cycle = 10;
nframes = gait_cycles * frames_per_cycle;

t_solution = linspace(0, gait_period*gait_cycles, nframes);
g_solution = deval(solution, mod(t_solution, gait_period));

gait_displacement = vec_to_mat_SE2(deval(solution, gait_period));


trajectory = zeros(size(g_solution));

for frame = 1:nframes
    t = t_solution(frame);
    gaits_finished = floor(t/gait_period);
    
    accumulated_displacement = gait_displacement ^ gaits_finished;
    g = accumulated_displacement * vec_to_mat_SE2(g_solution(:, frame));
    
    
    trajectory(:, frame) = mat_to_vec_SE2(g);
    
end
%figure()
plot(trajectory(1,:), trajectory(2, :))
axis equal

end

