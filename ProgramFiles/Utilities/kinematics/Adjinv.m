function Adjoint_inverse_action = Adjinv(g)
% Inverse adjoint action of the group element g

%invTgLginv = ThLg(g,[0;0;0]);

%TeRg =(TgRh([0;0;0],g));

Adjoint_inverse_action = TgLginv(g)*TeRg(g);

end