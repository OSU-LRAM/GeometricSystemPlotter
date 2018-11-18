function Adjoint_action = Adj(g)
% Adjoint action of the group element g

% Convert from matrix representation to column if needed.
if numel(g) == 9
    g = mat_to_vec_SE2(g);
end

Adjoint_action = TgRginv(g)*TeLg(g);%invTgRginv\TeLg;

if isa(Adjoint_action,'sym')
    Adjoint_action = simplify(Adjoint_action,'steps',50);
end

end