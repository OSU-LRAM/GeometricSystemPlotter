function valid = validate_dMalpha_dalpha_calc()
% Define jacobian using values from Ross's book for rotary-rotary arm
syms l1 l2 a1 a2 m1 m2 I1 I2 da1 da2 real;
J1 = [0 0; l1/2 0; 1 0];
J2 = [l1*sin(a2) 0; l2/2+l1*cos(a2) l2/2; 1 1];
J = {J1, J2};

% Derivative of J with respect to the shape variables
dJ1da1 = sym([0 0; 0 0; 0 0]);
dJ1da2 = sym([0 0; 0 0; 0 0]);
dJ2da1 = sym([0 0; 0 0; 0 0]);
dJ2da2 = sym([l1*cos(a2) 0; -l1*sin(a2) 0; 0 0]);
% Check the partials against the jacobian derivative calculator
dJdq = fixed_jacobian_derivative(J);
test1 = isequaln(dJ1da1,cell2sym(dJdq{1}(:,1)'));
test2 = isequaln(dJ1da2,cell2sym(dJdq{1}(:,2)'));
test3 = isequaln(dJ2da1,cell2sym(dJdq{2}(:,1)'));
test4 = isequaln(dJ2da2,cell2sym(dJdq{2}(:,2)'));
if ~(test1 && test2 && test3 && test4)
    error('Mismatch in Jacobian derivative calculation')
end

% Local inertia tensors
u1 = diag([m1 m1 I1]);
u2 = diag([m2 m2 I2]);

% % Check the mass matrix - debugging only
% M = J1'*u1*J1 + J2'*u2*J2;
% simplify(M)

% Find what the partial_mass_matrix function says are the partial mass
% terms
dMda = partial_mass_matrix(J,dJdq,{u1, u2},[]);
% Check these against the book's asserted values
dMda1 = sym([0 0; 0 0]);
dMda2 = -m2*l1*l2*sin(a2)/2*[2 1; 1 0];
test1 = isequaln(dMda{1},dMda1);
test2 = isequaln(dMda{2},dMda2);
if ~(test1 && test2)
%     error('Mismatch in dMdalpha calculation')
    valid = false;
else
    valid = true;
end

shape = [a1, a2];
dshape = [da1, da2];
C = calc_coriolis_matrix(dMda,shape,dshape);
% Check to see if the output of calc_coriolis_matrix is the same as the
% book's
Cbook = -m2*l1*l2*sin(a2)/2*[2*da1*da2+da2^2;da1*da2] + ...
    m2*l1*l2*sin(a2)/2*[0; da1^2 + da1*da2];
Cbook = simplify(Cbook);
if ~isequaln(C,Cbook)
    valid = false;
else
    valid = true;
end
end