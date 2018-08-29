function LA = Tg1Lg2(g,~)

theta = g(3);

LA = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];

end