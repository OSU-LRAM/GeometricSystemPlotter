function [h, J, J_full] = backbone_from_curvature_bases(kappa_basis_input,r,L)
% Integrate backbone geometry from a set of curvature basis functions and
% weighting factors on those basis functions. The backbone is
% of total length L, placed such that its midpoint is at the origin and its
% midpoint tangent is aligned with the x axis.
% 
%
% Inputs: 
%
% kappa_basis_input: A single-column cell array of function handles, each
%   of which maps from the domain s = [-.5 .5] (normalized unit-length
%   backbone) to the curvature (dtheta/ds) of the backbone at this point.
%
%   Curvature basis functions should be vectorized (so that if given a row
%   vector of s values, they return a row vector of curvatures.
%
%   Faster performance is achieved if the curvature basis functions return
%   a second *row* in their output, which contains the integral of the
%   curvature from the midpoint to that point on the body (the orientation
%   of the tangent line at that point on the body, relative to the tangent
%   line at the midpoint). If this integral is provided, best performance
%   will correspond to having it equal to zero when s is zero.
%
% r: A column of scalar values, which are used to scale the outputs of
%   their respective kappa functions, such that the total curvature of the
%   backbone is the r-weighted sum of the kappa_basis_input functions,
%
%      kappa(s) = sum( kappa_basis{i}(s) * r(i) )
%
% L: The actual length of the backbone, in "real world" units. This value
%   scales the geometry of the integrated backbone shape
%
% Outputs:
%
% h is a function on the domain s = [-.5 .5] (normalized unit-length
%   backbone) that returns the position and orientation (in world-units) of
%   the frame on the backbone that is at point s along the backbone (in
%   backbone units)
%
% J is the Jacobian of h, so that its product with a shape velocity gives
%   the rate of change of h -- the "world velocity" of the backbone point
%   taking the midpoint-tangent frame as "fixed"
%
% J_full is the Jacobian mapping from system-body-velocity and shape
%   velocity to the body velocity (relative to a non-moving frame) of the
%   point on the body

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%Input checking

%make sure kappa_basis is a column of function handles
if size(kappa_basis_input,2) ~= 1 || numel(size(kappa_basis_input)) ~= 2
	if ~isa(kappa_basis_input,'cell') && all(cellfun('isclass',kappa_basis_input,'function_handle'))
		error('kappa_basis_input should be a cell column vector of function handles')
	else
		error('kappa_basis_input should be a column vector')
	end
end

%make sure r is a vector of scalars
if min(size(r)) ~= 1 || numel(size(r)) ~= 2
	if ~isa(r,'numeric')
		error('r should be a numeric vector')
	else
		error('r should be a vector')
	end
end
r = r(:);   % Make sure this is a column
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

% Set the integration limits. Internally, everything uses [-0.5 0.5], and
% then scales by the length factor at the end.
all_limits = [-0.5 0.5];

%%%%%%%%
% If kappa_basis does not already include an integral term, add it, and
% normalize to the body length

% Make a cell array for the actual kappa_basis used by the rest of the code
kappa_basis = cell(size(kappa_basis_input,1),1);

for i = 1:length(r)
	
	% check if the curvature basis provides an integral function
	
	switch length(kappa_basis_input{i}(0))
		
		% no integral provided, use ODE45
		case 1
			
			%sol = ode45(@(s,theta) kappa_basis_input{i}(s),int_limit,0);
			%kappa_basis{i} = @(s) [kappa_basis_input{i}(s); deval(sol,s)-deval(sol,0)]; % Make zero-tangent at center
			
			% Piecewise integration
			sol = ode_multistart(@ode45,@(s,theta) kappa_basis_input{i}(s),all_limits,0,0);
			kappa_basis{i} = @(s) [kappa_basis_input{i}(s); sol(s)];
				
		% integral provided
		case 2
			
			kappa_basis{i} = kappa_basis_input{i};
			
		otherwise
			
			error('Basis gives unexpected output')
		
	end
	
end


%%%%%%
% Several of the anonymous functions below are admittedly a mess -- they're
% what I needed in order to get vectorized behavior from the output
% functions

%%%%%%%%%%%%%%%%%
% Get the xytheta coordinates of points on the body -- theta coordinates are
% handled separately, to make it easier to handle
% delta-functions and the like

%summed integral of the curvature functions
theta_fun = @(s) sum(cell2mat(cellfun( @(k,r) theta_fun_helper(k,r,s), kappa_basis,num2cell(r),'UniformOutput',false)),1);


% SE2 integrate along the backbone with known theta
%locus_sol = ode45(@(s,h) locus_helper(s,theta_fun),int_limit,[0;0]);
locus_sol = ode_multistart(@ode45, @(s,h) locus_helper(s,theta_fun),all_limits,0,[0;0]);






% Return the locus data as 3xn matrix with x,y,theta at each of n points
%h = @(s) [L*(deval(locus_sol,torow(s)/L)-deval(locus_sol,0*torow(s))); theta_fun(torow(s)/L)];
h_norm = @(s) [locus_sol(torow(s)); theta_fun(torow(s))]; % in normalized internal coordinates
h = @(s) [L*locus_sol(torow(s)); theta_fun(torow(s))]; % with actual length

%%%%%%%%%%%%%%
% Get the jacobian to body point velocities

if nargout > 1

    % Calculate the theta component separately once again, to accomodated delta
    % functions

    % Summed integral of the basis functions
    J_theta_fun = @(s) cell2mat(permute(cellfun( @(k) toz(theta_fun_helper(k,1,s)), kappa_basis,'UniformOutput',false),[2 1]));

    % SE2 integration for velocities
    %jacobian_sol = ode45(@(s,J) J_helper(s,J,kappa_basis,h(s)),L*int_limit,zeros(2*length(kappa_basis),1));
    jacobian_sol = ode_multistart(@ode45,@(s,J) J_helper(s,J,kappa_basis,h_norm(s)),all_limits,0,zeros(2*length(kappa_basis),1));

    % Concatenate xy and theta jacobians.
    J = @(s) cat(1,reshape(L*jacobian_sol(toz(s)),2,[],length(s)),J_theta_fun(s));

end

%%%%%%%
% Get the full Jacobian from system body velocity and shape velocity to
% body velocity of the point on the body, (relative to fixed world frame)
if nargout > 2
    
    J_full = @(s) [Adjinv(h(s)) TgLginv(h(s))*J(s)];
    
end
    
end

% Get the derivative of the locus along the backbone length
function dh = locus_helper(s,theta_fun)

	% Build the rotation matrix
	theta = theta_fun(s);
	R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
		
	% Rotate the basic tangent vector to point along the backbone
	dh = R * [1; 0];
	
end

% get the derivative of the jacobian along the backbone length
function dJ = J_helper(s,J,kappa_basis,h) %#ok<INUSL>

	dJ = zeros(2,length(kappa_basis));
	
	for i = 1:length(kappa_basis)
		
		int_basis = theta_fun_helper(kappa_basis{i},1,s);
		
		dJ(:,i) = [-sin(h(3)) -cos(h(3));
			cos(h(3)) -sin(h(3))] * [int_basis; 0];%* [[0 1]*(kappa_basis{i}(s)-kappa_basis{i}(0)); 0];
		
	end
	
	%Reshape into a vector for ODE
	dJ = dJ(:);

end





function angle = theta_fun_helper(k,r,s) 

	angle_basis = (k(torow(s)));
	angle = r*angle_basis(2,:);
	
end