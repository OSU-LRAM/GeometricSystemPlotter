function [h, J] = backbone_from_curvature_bases(kappa_basis_input,r,L)
% Build up the backbone locus and velocity one-form at a given shape
% defined by the curvature kappa = sum(r(i)*kappa_basis(i). The backbone is
% of total length L, extending in each direction from the point (0,0) with
% orientation 0 (i.e. along the x axis). Curvature basis functions kappa
% should be defined with respect to a total body length of 1, should be
% vectorized, and can provide an explicit integral of their function as a
% second row of the output. If this integral is provided, best performance
% will correspond to having it pass through (0,0)


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%Input checking -- this probably isn't really necessary, but I want to
%practice it

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
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

% Set the integration limits. Internally, everything uses [-0.5 0.5], and
% then scales by the length factor at the end.
int_limit = [-0.5 0.5];


%%%%%%%%%%%%%%%%
% Sort the discontinuities into the integration limit vector
all_limits = sort([int_limit 0 discont]);


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
h = @(s) [L*locus_sol(torow(s)/L); theta_fun(torow(s)/L)]; % with actual length

%%%%%%%%%%%%%%
% Get the jacobian to body point velocities

if nargout == 2

    % Calculate the theta component separately once again, to accomodated delta
    % functions

    % Summed integral of the basis functions
    J_theta_fun = @(s) cell2mat(permute(cellfun( @(k) toz(theta_fun_helper(k,1,s)), kappa_basis,'UniformOutput',false),[2 1]));

    % SE2 integration for velocities
    %jacobian_sol = ode45(@(s,J) J_helper(s,J,kappa_basis,h(s)),L*int_limit,zeros(2*length(kappa_basis),1));
    jacobian_sol = ode_multistart(@ode45,@(s,J) J_helper(s,J,kappa_basis,h_norm(s)),all_limits,0,zeros(2*length(kappa_basis),1));

    % Concatenate xy and theta jacobians.
    J = @(s) cat(1,reshape(L*jacobian_sol(toz(s/L)),2,[],length(s)),J_theta_fun(s/L));

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


function X = torow( X )

	X = X(:).';
    
end

function Y = toz( X )

	%first, make X a row
	X = torow(X);
	
	% then, make it a 3rd-dimension column
	Y(1,1,:) = X;
    
end

function angle = theta_fun_helper(k,r,s) 

	angle_basis = (k(torow(s)));
	angle = r*angle_basis(2,:);
	
end