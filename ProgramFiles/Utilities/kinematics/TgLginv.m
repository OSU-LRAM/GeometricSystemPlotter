function LA = TgLginv(g)
% Left lifted action of g bringing vectors to the origin/Lie algebra

% Prevent Matlab from playing tricks with imaginary numbers on symbolic
% inputs and from complaining about assumptions on constants
if isa(g,'sym')
    assume(symvar(g),'real');
    warning('off','symbolic:sym:sym:AssumptionsOnConstantsIgnored')
end


% Convert from matrix representation to column if needed.
if numel(g) == 9
    g = mat_to_vec_SE2(g);
end


theta = g(3);

LA = [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];

end