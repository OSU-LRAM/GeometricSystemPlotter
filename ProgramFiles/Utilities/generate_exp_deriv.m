% computes the time-derivative of the exponential map
% ouputs:
    % f: anonymous function, accepting (g, d/dt g) as vectors
function f = generate_exp_deriv
    % define symvars for beta, derivatives
    beta_vec = sym('beta', [1 3]);
    d_beta_vec = sym('d_beta', [1 3]);
    assume(beta_vec, 'real');
    assume(d_beta_vec, 'real');
    % generate matrices
    beta_mat = vec_to_mat_SO3(beta_vec);
    d_beta_mat = vec_to_mat_SO3(d_beta_vec);
    % compute integral
    t = sym('t');
    grand = expm(-t*beta_mat) * d_beta_mat * expm(t*beta_mat);
    p_exp_mat = real(int(grand, t, 0, 1)); %ignore zero-valued complex part
    f = matlabFunction(p_exp_mat,...
                   'file', 'exp_deriv',...
                   'vars', {beta_vec, d_beta_vec});
end