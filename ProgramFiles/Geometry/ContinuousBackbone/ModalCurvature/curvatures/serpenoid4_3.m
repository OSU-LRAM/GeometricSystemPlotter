function kappa = serpenoid4_3(s,omega)
%sinusoidal curvature

if ~exist('omega','var')
    omega = 1;
end

kappa = omega*cos(4*pi*omega*s);