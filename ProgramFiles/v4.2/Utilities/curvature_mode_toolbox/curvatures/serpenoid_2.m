function kappa = serpenoid_2(s,omega)
%sinusoidal curvature

if ~exist('omega','var')
    omega = 1;
end

kappa = sin(2*pi*omega*s);