function LA = TgLginv(g)
% Left lifted action of g bringing vectors to the origin/Lie algebra

% Convert from matrix representation to column if needed.
if numel(g) == 9
    g = mat_to_vec_SE2(g);
end


theta = g(3);

LA = [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];

end