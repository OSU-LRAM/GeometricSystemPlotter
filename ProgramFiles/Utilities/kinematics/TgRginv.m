function RA = TgRginv(g)
% Right lifted action of g bringing vectors to the origin/Lie algebra

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


xdelta = g(1);
ydelta = g(2);

RA = [1 0 ydelta;
    0 1 -xdelta;
    0 0 1];

end