function gait_function = generate_1D_gait(amplitude, center, direction)
%GENERATE_1D_GAIT Summary of this function goes here
%   Detailed explanation goes here

% entry 1: anonymous function for a(t)
% entry 2: anonymous function for adot(t) (derivative of a(t))

gait_function = {@(t) center + amplitude * sin(t * direction), ...
                 @(t) direction * amplitude * cos(t * direction)};

end

