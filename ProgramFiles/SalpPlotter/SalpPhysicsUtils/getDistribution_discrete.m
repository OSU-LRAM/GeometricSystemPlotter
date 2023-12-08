function B = getDistribution_discrete(J_full)
% Get distribution matrix for a set of forces
%
% Inputs:
%
%   J_full: Cell array with full Jacobians from body velocity of the base
%           frame and joint angle velocities to body velocity of each
%           thruster location
%
% Outputs:
%
%   B: Distribution matrix for a given set of Jacobians

J_dual_full = cellfun(@transpose, J_full, 'UniformOutput', false);
B = [J_dual_full{:}];

