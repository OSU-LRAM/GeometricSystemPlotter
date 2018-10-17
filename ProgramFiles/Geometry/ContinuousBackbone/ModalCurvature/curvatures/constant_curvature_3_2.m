function kappa = constant_curvature_3_2(s)
%Curvature for first half of the system

kappa = zeros(size(s));
% kappa(s>-1/6 && s<1/6) = 1;
if s>-1/6
    if s<1/6
        kappa=1;
    end
end
       