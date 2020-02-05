function [t_solution, g_solution] = asym_solve_gait(s, gait,system_map)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t0 = 0;
tf = 2*pi;
g0 = [0 0 0];
[t_solution, g_solution] = ode45(@(t, g) calc_gdot(t, g, s, system_map, gait), [t0, tf], g0);

end

function gdot = calc_gdot(t, g, s, system_map, gait)
    gdot = TeLg(g) * apply_piecewise_system(s, system_map, gait{1}(t), gait{2}(t));
end