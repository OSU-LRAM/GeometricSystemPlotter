function [body_velocity] = apply_piecewise_system(s, piecewise_system_map, shape, shapechange)
%APPLY_PIECEWISE_SYSTEM Summary of this function goes here
%   Detailed explanation goes here

% For now, determine dir according to the sign of a and adot
% using results found for two-link.

% pick which system to use for this a/adot
A = get_correct_A(s, shape, shapechange);

% gcirc right
body_velocity = -A * shapechange;

end