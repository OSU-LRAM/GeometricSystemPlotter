function RA = Tg1Rg2(g1,g2)

theta0 = g1(3);
xdelta = g2(1);
ydelta = g2(2);

RA = [1 0 -(xdelta*sin(theta0) + ydelta*cos(theta0));
    0 1 (xdelta*cos(theta0) - ydelta*sin(theta0));
    0 0 1];

end