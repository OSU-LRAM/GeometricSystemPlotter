function kappa = constant_curvature_3_3(s)
%Curvature for first half of the system

kappa = zeros(size(s));
kappa(s>1/6) = 1;