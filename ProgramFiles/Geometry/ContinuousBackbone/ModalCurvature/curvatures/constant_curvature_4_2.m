function kappa = constant_curvature_4_2(s)
%Curvature for first half of the system

kappa = zeros(size(s));
if s>-1/4
    if s<0
        kappa=1;
    end
end