function [A] = get_A_of_alpha(s, grid)
%GET_A_OF_ALPHA Summary of this function goes here
%   Detailed explanation goes here

% connection_grid = [a_grid index, xyt, increasing or decreasing adot]
A_positive = zeros([length(grid),3]);
A_negative = zeros([length(grid),3]);
for i=1:length(grid)
    A_positive(i,:) = -get_correct_A(s, grid(i), 1);
    A_negative(i,:) = -get_correct_A(s, grid(i), -1);
end

A_difference = A_positive - A_negative;

A.positive = A_positive;
A.negative = A_negative;
A.difference = A_difference;

end

