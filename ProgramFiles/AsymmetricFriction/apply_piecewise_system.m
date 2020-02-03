function [body_velocity] = apply_piecewise_system(s, piecewise_system_map, shape, shapechange)
%APPLY_PIECEWISE_SYSTEM Summary of this function goes here
%   Detailed explanation goes here

% For now, determine dir according to the sign of a and adot
% using results found for two-link.

% pick which system to use for this a/adot
opening = sign(shape) == sign(shapechange);
if opening
    dir = 3 - 1;
else
    dir = 2 - 1;
end
% dir = piecewise_system_map(i, j) - 1;

backwards = [mod(dir,2), floor(dir/2)]; % this is a hacky way to do all 4 direction combos with binary

% get A(alpha)
[A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry, s.physics, shape, backwards);

% gcirc right
body_velocity = A * shapechange;

end

