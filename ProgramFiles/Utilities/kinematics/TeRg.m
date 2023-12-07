function RA = TeRg(g)
% Right lifted action of g acting on vectors at the origin/in the Lie
% algebra

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

RA = vec_to_mat_SE2([-ydelta,xdelta,0]);

end