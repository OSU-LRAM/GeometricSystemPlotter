function Adjoint_inverse_action = Adjinv(g)
% Inverse adjoint action of the group element g

% Convert from matrix representation to column if needed.
if numel(g) == 9
    g = mat_to_vec_SE2(g);
end

Adjoint_inverse_action = TgLginv(g)*TeRg(g);

end