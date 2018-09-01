function RA = TgRginv(g)
% Right lifted action of g bringing vectors to the origin/Lie algebra

% Convert from matrix representation to column if needed.
if numel(g) == 9
    g = mat_to_vec_SE2(g);
end


xdelta = g(1);
ydelta = g(2);

RA = [1 0 ydelta;
    0 1 -xdelta;
    0 0 1];

end