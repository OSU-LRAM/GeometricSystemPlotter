function A = get_correct_A(s, shape, shapechange)
%GET_CORRECT_A Summary of this function goes here
%   Detailed explanation goes here

opening = sign(shape) == sign(shapechange);

if opening
    dir = 3 - 1;
else
    dir = 2 - 1;
end

backwards = [mod(dir,2), floor(dir/2)]; % this is a hacky way to do all 4 direction combos with binary

% get A(alpha)
[A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry, s.physics, shape, backwards);
end

