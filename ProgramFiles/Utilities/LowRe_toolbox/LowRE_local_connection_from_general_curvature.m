function [A, h, J,Omega] = LowRE_local_connection_from_general_curvature(curvdef,cparams,L,c,drag_ratio)
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

% First, get the backbone locus and Jacobian functions
[h, J] = backbone_from_general_curvature(curvdef,cparams,L);

% Now integrate to get the pfaffian
%Omega_sol = ode45( @(s,Omega) connection_helper(s,h(s),J(s)),int_limit,zeros(3,3+length(r)));
Omega_sol = ode_multistart(@ode45,@(s,Omega) connection_helper(s,h(s),J(s),c,drag_ratio),int_limit,int_limit(1),zeros(3,3+length(cparams)));

Omega = reshape(Omega_sol(int_limit(end)),3,[]);

A = Omega(:,1:3)\Omega(:,4:end);

end

function dOmega = connection_helper(s,h,J,c,drag_ratio) %#ok<INUSL>
% Calculate the derivative of the local connection as it's built up along
% the backbone

	% Midpoint-tangent-frame force jacobian
	gdot_to_xi_local = [cos(h(3)) sin(h(3)) 0;
			-sin(h(3)) cos(h(3)) 0;
			0 0 1];
		
	% Local drag, based on unit longitudinal drag, lateral according to the ratio, no local
	% torsional drag, multiplied by drag coefficient
	xi_local_to_F_local = [-1 0 0;0 -drag_ratio 0;0 0 0]*c;
	
	% Force local element applies at midpoint-tangent-frame
	F_local_to_F_midpoint = [cos(h(3)) -sin(h(3)) 0;
			sin(h(3)) cos(h(3)) 0;
			0 0 1];
		
	% Moment around midpoint frame induced by force
	F_midpoint_to_FM_midpoint = [1 0 0; 0 1 0; -h(2) h(1) 1];
	
	% shape component of pfaffian
	omega2 = F_midpoint_to_FM_midpoint * F_local_to_F_midpoint...
		* xi_local_to_F_local * gdot_to_xi_local * J;
	
		% body velocity component of pfaffian
	omega1 = F_midpoint_to_FM_midpoint * F_local_to_F_midpoint...
		* xi_local_to_F_local * gdot_to_xi_local * [1 0 -h(2); 0 1 h(1); 0 0 1];
	
	dOmega = [omega1 omega2];
	dOmega = dOmega(:);
	
% 	% local contribution to the local connection
% 	dA = omega1\omega2;
% 	
% 	% reshape into a vector
% 	dA = dA(:);

	

end

function X = torow( X )
%
% TOROW Converts a vector or a matrix into a row vector.
%   If input is already a row vector, it is returned with no change.
%   If input is a column vector, it is converted into a row vector and
%   returned.
%   If input is a matrix, each column is converted into a row, and all
%   resulting rows are placed in series into a single row which is
%   returned.
%
% Input:
%   X - input vector or matrix
%
% Output:
%   X - row vector
%
% Examples:
%   torow([ 0 1 2 3 ])
%       returns [ 0 1 2 3 ]
%   torow([ 0 1 2 3 ]')
%       returns [ 0 1 2 3 ]
%   torow([ 0 1; 2 3 ])
%       returns [ 0 2 1 3 ]
%   torow([ 0 1; 2 3 ]')
%       returns [ 0 1 2 3 ]
%
% Author:	Tashi Ravach
% Version:	1.0
% Date:     07/07/2010
%

    % check if input is a vector
    [ m, n ] = size(X);
    if m == 1
        return % input is already a row vector with n columns
    elseif n==1
        X = X'; % input is converted from column vector to row vector
    elseif (m*n>n) || (m*n>m)
        X = X(:)'; % input is converted from matrix to row vector by column
    else
        X = []; % input is unknown and an empty output is returned
    end
    
end

function Y = toz( X )

	%first, make X a row
	X = torow(X);
	
	% then, make it a 3rd-dimension column
	Y(1,1,:) = X;
    
end