function [kappa, joint_loc] = discrete_3_joint_3(s)
% Curvature and integral of curvature for a discrete joint 3/4 of the
% way down the body, specified relative to a total body length of 1 unit,
% centered at 1/2
%
% By returning joint_loc as a second output, this basis function provides
% an integral of the curvature that is more robust and faster to compute
% than directly integrating a Dirac delta function.

%joint is at 1/2 way from midpoint to tail end
joint_loc = 1/4;

% Curvature 
kappa = [dirac(s-joint_loc); heaviside(s-joint_loc)-heaviside(0-joint_loc)];