function solution = asym_solve_gait(s, gait)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
t0 = 0;
tf = 2*pi;
g0 = [0 0 0];
solution = ode45(@(t, g) calc_gdot(t, g, s, gait), [t0, tf], g0);
end

function gdot = calc_gdot(t, g, s, gait)
    gdot = TeLg(g) * apply_piecewise_system(s, gait{1}(t), gait{2}(t));
end