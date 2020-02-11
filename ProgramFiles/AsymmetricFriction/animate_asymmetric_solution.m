function F = animate_asymmetric_solution(system, solution, gait)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

gait_period = 2*pi;
gait_cycles = 3;
frames_per_cycle = 25;
nframes = gait_cycles * frames_per_cycle;

t_solution = softspace(0,gait_period * gait_cycles, nframes);
g_solution = deval(solution, mod(t_solution, gait_period));
gait_displacement = deval(solution, gait_period);

% record movie
figure()
axis manual
ax = gca;
ax.NextPlot = 'replaceChildren';
xlim([-1, 1]);
ylim([-1, 1]);

F = struct('cdata',[],'colormap',[]); % reset length of struct
F(nframes) = struct('cdata',[],'colormap',[]);
for frame = 1:nframes
    t = t_solution(frame);
    gaits_finished = floor(t/gait_period);
    
    
    % bug: this should be a matrix transformation, not an addition
    accumulated_displacement = gait_displacement * gaits_finished;
    g = g_solution(:, frame) + accumulated_displacement;
    
    
    B = vec_to_mat_SE2(g) * fat_chain(system.geometry, gait{1}(t));
    cla
    hold on
    patch(B(1,:),B(2,:),[1 0 0])
    drawnow
    F(frame) = getframe;
end

end

