function [e2_mean, tissot, meantau] = mean_distortion_lite(metric,dist_function,x1,x2,y1,y2,Jac)
% Calculate the mean distortion of the input metric for the specified
% distortion function. x and y inputs are as quad2d, and J is the jacobian
% matrix of which the determinant is the area scaling.
%
% The second output is the Tissot distortion
%
% This "lite" version replaces the quadrature calls with simple array
% summing

	% If there is no jacobian specified, make it unit
	if ~exist('Jac','var')
		Jac = eye(2);
	end

	% If the jacobian is not a function handle, make a dummy function
	if ~isa(Jac,'function_handle')
		Jac = @(x,y) Jac;
	end
	
	% Make a grid
	res = 20;	
	[xgrid,ygrid] = meshgrid(linspace(x1,x2,res), linspace(y1,y2,res));
	
	%Find the minimum stretch in the domain
	mintau = find_min_stretch(metric,x1,x2,y1,y2);
	
	%find the maximum stretch in the domain
	maxtau = find_max_stretch(metric,x1,x2,y1,y2);

	% Calculate the area of the domain
	%A = arrayquad2d(@(x,y) det(Jac(x,y)),x1,x2,y1,y2);
	dA = arrayfun(@(x,y) det(Jac(x,y)),xgrid,ygrid);
	A = sum(dA(:));
	
    %Get the geometric mean principle stretch in the domain
    dtau_tot = arrayfun(@(x,y) geomean(tissot_stretches(metric,x,y))...
		...*det(Jac(x,y))
        ,xgrid,ygrid);
	%tau_tot = sum(dtau_tot(:));


	% Find the geometric mean stretch as the geometric mean of the
	% local geometric means weighted by the area they represent
	meantau = exp((dA(:)'*log(dtau_tot(:)))/A);
	
	% Calculate the integrated distortion
% 	e2_tot = arrayquad2d(@(x,y) metric_distortion(metric,dist_function,x,y,meantau)...
% 		*det(Jac(x,y)),x1,x2,y1,y2);
	de2_tot = arrayfun(@(x,y) metric_distortion(metric,dist_function,x,y,meantau)...
		*det(Jac(x,y)),xgrid,ygrid);
	e2_tot = sum(de2_tot(:));
	


	% Average distortion
	e2_mean = e2_tot/A;
	
	% Maximum distortion ratio
	tissot = maxtau/mintau;

end

function output = arrayquad2d(fun,x1,x2,y1,y2)
% Wrapper for quad2d that evaluates the integrand by way of arrayfun

	output = quad2d(@(X,Y) arrayfun(@(x,y) fun(x,y),X,Y),x1,x2,y1,y2);
	
end

function tau = tissot_stretches(M,x,y)
% Find the principle stretches of the metric at the given point as a cell
% array

	[u,s,v] = svd(inv(M(x,y))); %#ok<NASGU,ASGLU>
	
	tau = diag(sqrt(s));
	
end

function e2 = metric_distortion(metric,dist_function,x,y,meantau)
% Calculate the distortion of the metric at a specified point, with
% stretches specified relative to the mean stretch meantau

	tau = tissot_stretches(metric,x,y)/meantau;
	e2 = dist_function(tau(1),tau(2));
	
end
	
function mintau = find_min_stretch(metric,x1,x2,y1,y2)
% Find the minimum stretch in the metric

	objective_func = @(X) find_min_stretch_helper(metric,X(1),X(2));
	
	mintau = fmincon_quadbound4(objective_func,x1,x2,y1,y2);
	
end

function maxtau = find_max_stretch(metric,x1,x2,y1,y2)
% Find the maximum stretch in the metric

	objective_func = @(X) - find_max_stretch_helper(metric,X(1),X(2));
	
	maxtau = - fmincon_quadbound4(objective_func,x1,x2,y1,y2);
	
end

function mintau = fmincon_quadbound4(objective_func,x1,x2,y1,y2)
% Get the find the minimum value of a 2d function over a region with bounds
% specified as in quad2d. 4 starting points are selected in an even
% distribution across the interior of the region

	% Convert y1 and y2 limits into functions if necessary
	if ~isa(y1,'function_handle')
		y1 = @(x) y1;
	end
	
	if ~isa(y2,'function_handle')
		y2 = @(x) y2;
	end
	
	% Build constraint definition for fmincon
	cA = [-1 0;1 0];
	cB = [-x1;x2];
	cAeq = [];
	cBeq = [];
	LB = [];
	UB = [];
	cNONLCON = @(X) yconst(X,y1,y2);
	cOptions = optimset('algorithm','interior-point','display','off');
	
	
	% Pick 9 starting points distributed evenly within x, and between the y
	% bounds for each x
	
	x_start_raw = linspace_interior(x1,x2,2)';
	x_start = repmat(x_start_raw,2,1);
	y_start = zeros(4,1);
	for i = 1:2
		y_start((1:2)+mod(i-1,3)*2) = linspace_interior(y1(x_start(i)),y2(x_start(i)),2);
	end
	
	% Find the minimum stretch from the various starting points
	mintau_raw = zeros(4,1);
	for i = 1:4
		

		
		[junk,mintau_raw(i)] = fmincon(objective_func,[x_start(i);y_start(i)],...
			cA,cB,cAeq,cBeq,LB,UB,cNONLCON,cOptions); %#ok<ASGLU>
		
	end
	
	% Get the minimum stretch over the different starting points
	mintau = min(mintau_raw);

end

function tau_1 = find_max_stretch_helper(M,x,y)
% Extract the smaller stretch

	tau = tissot_stretches(M,x,y);
	tau_1 = tau(1);

end

function tau_2 = find_min_stretch_helper(M,x,y)
% Extract the smaller stretch

	tau = tissot_stretches(M,x,y);
	tau_2 = tau(2);

end

function [C,Ceq] = yconst(X,y1,y2)
% Nonlinear constraints for fmincon that uses the y bound style from quad2d

	C = [y1(X(1)) - X(2);X(2) - y2(X(1))];
	
	Ceq = [];
	
end
