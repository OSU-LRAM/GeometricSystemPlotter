function kappa = serpenoid4_4(s,omega)
%sinusoidal curvature

if ~exist('omega','var')
    omega = 1;
end

kappa = omega*sin(4*pi*omega*s);