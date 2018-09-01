function LA = TeLg(g)
% Left lifted action of g acting on vectors at the origin/in the Lie
% algebra

% Convert from matrix representation to column if needed.
if numel(g) == 9
    g = mat_to_vec_SE2(g);
end

theta = g(3);

LA = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];

end