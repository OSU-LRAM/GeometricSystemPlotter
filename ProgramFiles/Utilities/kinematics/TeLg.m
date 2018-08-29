function LA = TeLg(g)
% Left lifted action of g acting on vectors at the origin/in the Lie
% algebra

theta = g(3);

LA = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];

end