function kappa = constant_curvature_4_1(s)
%Curvature for first half of the system

kappa = zeros(size(s));
kappa(s<-1/4) = 1;