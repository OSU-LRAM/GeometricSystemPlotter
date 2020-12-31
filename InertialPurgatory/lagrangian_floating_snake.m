clc
clear
syms m I alpha1 alpha2 dalpha1 dalpha2 c1 c2 s1 s2

M = [m 0 0;0 m 0; 0 0 I];   % Inertia Matrix
l = sym('l',[1,3]);         % Links length

dxb = sym('dxb',[1,3]);

g_circ = sym('g_circ',[1,3]);

% T = @(tt) [cos(tt) sin(tt) 0; -sin(tt) cos(tt) 0; 0 0 1];

Left_lifted_action = @(theta) [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];

Right_lifted_action = @(x,y) [1 0 -y; 0 1 x; 0 0 1];


% LA1 = [c1 -s1 0;s1 c1 0;0 0 1];
% LA2 = [c2 s2 0;-s2 c2 0;0 0 1];


% Find the body velocity for left link
f2 = Right_lifted_action(-l(2)/2,0)*transpose(g_circ);
h1 = Left_lifted_action(-alpha1)*(f2 + [0;0;-dalpha1]);
% h1 = LA1*(f2 + [0;0;-dalpha1]);
g1_circ = Right_lifted_action(-l(1)/2,0)*h1;

% Find the body velocity for left link
h2 = Right_lifted_action(l(2)/2,0)*transpose(g_circ);
f3 = Left_lifted_action(alpha2)*(h2 + [0;0;dalpha2]);
% f3 = LA2*(h2 + [0;0;dalpha2]);
g3_circ = Right_lifted_action(l(3)/2,0)*f3;

% Find the body velocity for middle link
g2_circ = transpose(g_circ);


var = [g_circ dalpha1 dalpha2];
var = [g_circ dalpha1 dalpha2;1:length(var)];

for j = 1:size(g1_circ,1)
    for i = 1:length(var)

        [c1,t1] = coeffs(g1_circ(j),[var(1,i)]);

        if (t1(1)) ~= 1
            J1(j,i) = c1(1);
        else
            J1(j,i) = 0;
        end

        [c2,t2] = coeffs(g2_circ(j),[var(1,i)]);

        if t2(1) ~= 1
            J2(j,i) = c2(1);
        else
            J2(j,i) = 0;
        end

        [c3,t3] = coeffs(g3_circ(j),[var(1,i)]);

        if (t3(1)) ~= 1
            J3(j,i) = c3(1);
        else
            J3(j,i) = 0;
        end

    end
end

M_new = 0.5*transpose(J1)*M*J1 + 0.5*transpose(J2)*M*J2 +...
    0.5*transpose(J3)*M*J3;

II = M_new(1:3,1:3);
IA = M_new(1:3,4:5);

A = simplify(II\IA);

% Find the Metric
M = transpose([-A; eye(2)])*M_new*[-A;eye(2)];

% ----------------------------------------------------------------------- %

% % Lagrangian Equation
% L = 0.5*transpose(g1_circ)*M*g1_circ + 0.5*transpose(g2_circ)*M*g2_circ +...
%     0.5*transpose(g3_circ)*M*g3_circ;
% 
% % Simplify the equation
% L = simplify (L);
% 
% J1 = jacobian(L,[g_circ dalpha1 dalpha2]);
% 
% 
% J2 = jacobian(J1,[g_circ dalpha1 dalpha2]);
% 
% for i =1:length(J2)
%     J2(i,i) = 0.5*J2(i,i);
% end
% 
% II1 = J2(1:3,1:3);
% IA1 = J2(1:3,4:5);
% 
% A1 = simplify(II1\IA1)
% 
% var = [g_circ dalpha1 dalpha2];
% var = [g_circ dalpha1 dalpha2;1:length(var)];
% 
% % M = zeros(5,5);
% for i= 1:length(var)
% 
%     [c,t] = coeffs(L,[var(1,i)]);
%     
%     if length(t) > 2
%         M_k(i,i) = c(1);
%     
%         var_temp = var;
%         var_temp(:,i) = [];
%         for j=1:length(var_temp)
%             [cs,ts] = coeffs(c(2),[var_temp(1,j)]);
%             
%             if length(ts) > 1
%                 M_k(i,var_temp(2,j)) = cs(1);
%             else
%                 M_k(i,var_temp(2,j)) = 0;
%             end
% 
%         end
%         
%     end
%         
% end
% 
% 
% II = M_k(1:3,1:3);
% IA = M_k(1:3,4:5);
% 
% A = simplify(II\IA);
% 
% % Find the Metric
% M = transpose([-A; eye(2)])*M_k*[-A;eye(2)];

f1 = matlabFunction(A1, 'file', 'A_func');

f2 = matlabFunction(M, 'file', 'M_func');

f2 = matlabFunction(M, 'file', 'M_func');

