function kappa = constant_curvature_2(s)
%Curvature for second half of the system

kappa = zeros(size(s));
kappa(s>0) = 1;