function kappa = serpenoid4_2(s,omega)
%sinusoidal curvature

if ~exist('omega','var')
    omega = 1;
end

kappa = omega*sin(2*pi*omega*s);