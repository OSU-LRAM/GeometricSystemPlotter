function kappa = constant_curvature_1(s)
%Curvature for first half of the system

kappa = zeros(size(s));
kappa(s<0) = 1;