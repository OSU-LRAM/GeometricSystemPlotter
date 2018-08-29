function Ladji = adji(g)

invTgLginv = ThLg(g,[0;0;0]);

TeRg =(TgRh([0;0;0],g));

Ladji = invTgLginv\TeRg;

end