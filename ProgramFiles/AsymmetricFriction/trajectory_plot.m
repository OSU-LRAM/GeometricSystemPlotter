function trajectory_plot(solution)
%TRAJECTORY_PLOT Summary of this function goes here
%   Detailed explanation goes here

time = linspace(0, 2* pi, 20);
dif_eq = deval(solution, time);
dif_eq_trajectory = dif_eq(1:2, :);
figure()
plot(dif_eq_trajectory(1,:), dif_eq_trajectory(2, :))
axis equal

end

