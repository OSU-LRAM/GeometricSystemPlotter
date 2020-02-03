function [body_velocity] = apply_piecewise_system(s, piecewise_system_map, adot_grid, i, j)
%APPLY_PIECEWISE_SYSTEM Summary of this function goes here
%   Detailed explanation goes here

a_grid = s.grid.eval{1};

% pick which system to use for this a/adot
dir = piecewise_system_map(i, j) - 1;
backwards = [mod(dir,2), floor(dir/2)]; % this is a hacky way to do all 4 direction combos with binary

shape = a_grid(i);
shapechange = adot_grid(j);

% get A(alpha)
[A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry, s.physics, shape, backwards);

% gcirc right
body_velocity = A * shapechange;

end

