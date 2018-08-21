function kappa = serpenoid_1(s,omega)
%sinusoidal curvature

if ~exist('omega','var')
    omega = 1;
end

kappa = omega*cos(2*pi*omega*s);