function Adjoint_action = Adj(g)
% Adjoint action of the group element g

%TeLg = ThLg(g,[0;0;0]);

%invTgRginv =(TgRh([0;0;0],g));

Adjoint_action = TgRginv(g)*TeLg(g);%invTgRginv\TeLg;

end