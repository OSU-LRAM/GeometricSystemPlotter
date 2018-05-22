function [A, h, J] = HighRE_local_connection_from_curvature_bases(kappa_basis_input,r,L,density)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
% kappa_basis_input: cell vector of handles to basis functions
% r: coefficients for basis functions
% L: total length of swimmer
% c: drag per unit length
% drag_ratio: ratio of lateral to longitudinal drag

% Specified integration limits
int_limit = L*[-0.5 0.5];

if ~exist('density','var')
    density = @(s) 1;
end

% First, get the backbone locus and Jacobian functions
[h, J] = backbone_from_curvature_bases(kappa_basis_input,r,L);




% Now integrate to get the full kinetic energy tensor
M_full_sol = ode_multistart(@ode45,@(s,Mp) dM(s,Mp,h,J,density),int_limit,int_limit(1),zeros((3+length(r))^2,1));
M_full = reshape(M_full_sol(int_limit(end)),3+length(r),[]);

%Extract connection and metric from full kinetic energy tensor

Mg = M_full(1:3,1:3);
Mga = M_full(1:3,4:end);
Ma = M_full(3:end,4:end);

A = Mg\Mga;

%M = (A.'*Mg*A) - (A.'*Mga) - (Mga.'*A) + Ma;

end


function localJ = local_body_velocity_J(h,J)
% Calculate the Jacobian from shape parameter velocity to local tangential
% and normal velocity

	adjinv = [cos(h(3)) sin(h(3)) 0;
		-sin(h(3)) cos(h(3)) 0
        0 0 1]*...
        [1 0 -h(2)
        0 1 h(1)
        0 0 1];
    
    Linv = [cos(h(3)) sin(h(3)) 0;
		-sin(h(3)) cos(h(3)) 0
        0 0 1];
    
    localJ = [adjinv Linv*J];

end

function dMKE = dM(s,Mp,h,J,density) %#ok<INUSL>
% Calculate the local contribution to the power metric

	localJ = local_body_velocity_J(h(s),J(s));
	
% 	dMKE = localJ.'*density(s)*(diag([1,1,0]) + 0*diag([(0.25*s)^2,1*1*s^2,(1/8)*(s^2 - (0.25*s)^2)^2]))*localJ;
    phro = 1/pi;
    a = s;
    b = s/3;
    m = phro*pi*a*b;
    I = m*(a^2 + b^2)/4;
    m1f = phro*pi*b^2;
    m2f = phro*pi*a^2;
    If = (1/8)*phro*pi*(a^2 - b^2)^2;
    dMKE = localJ.'*density(s)*(diag([m,m,I]) + 0*diag([m1f, m2f, If]))*localJ;
	
	dMKE = dMKE(:);
	
end