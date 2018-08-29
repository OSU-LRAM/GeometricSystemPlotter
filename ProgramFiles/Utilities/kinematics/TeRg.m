function RA = TeRg(g)

xdelta = g(1);
ydelta = g(2);

RA = [1 0 -ydelta;
    0 1 xdelta;
    0 0 1];

end