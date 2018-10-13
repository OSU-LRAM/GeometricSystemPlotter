function kappa = constant_curvature_4_3(s)
%Curvature for first half of the system

kappa = zeros(size(s));

if s>0
    if s<1/4
        kappa=1;
    end
end