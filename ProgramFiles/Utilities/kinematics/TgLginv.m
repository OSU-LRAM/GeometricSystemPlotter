function LA = TgLginv(g)
% Left lifted action of g bringing vectors to the origin/Lie algebra

theta = g(3);

LA = [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];

end