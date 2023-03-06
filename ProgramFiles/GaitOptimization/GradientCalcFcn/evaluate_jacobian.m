function [jacobian,lineint,totalstroke,y,invert] = evaluate_jacobian(f,s,n,dimension,direction,lb,ub,~)
%%%%%%%%%%%%%
% This function calculates the gradient of displacement and cost function
% with respect to the coefficients obtained by 
% the fourier series parametrization
% 
% Inputs: 
% 
% y: Matrix containing the Fourier series coefficients
% s: System file which contains the connection vector field, CCF's and
%   metric data
% dimension: Indicates the number of shape variables of the system
% n: The number of points desired in a direct transcription parametrization
%   of the gaits
% direction: direction in which to optimize motion: 1-x, 2-y, 3-theta
% costfunction: costfunction type to optimize for
% 
% Outputs:
%
% jacobfourier: the struct which contains the gradient of the functions.
%   (jacobfourier.disp, jacobfourier.stroke, jacobfourier.eqi,
%   jacobfourier.total)
% lineint : the displacement for the gait
% totalstroke : the cost for the gait
%%%%%%%%%%%%%
if (isstruct(f))
    y = [f.phi_def{1}(linspace(0,T,n).'),f.phi_def{2}(linspace(0,T,n)).'];
    T = 2*pi;
else
    y = path_from_fourier(f,n,dimension);
    p = makeGait(f);
    T = 2*pi/f(end,1);
    y = y(1:end-1,:);
end
% Calculate displacement, cost and efficiency of a gait
% Note that, for inertial cost, cost is returned as the integral of torque
% squared, while for drag-based systems, cost is the path length of the
% gait
[~, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');
lineint=net_disp_opt(direction); % displacement produced in the chosen direction produced on executing the gait measured in the optimal coordinates 

% Assign value for totalstroke, i.e. the cost metric used for efficiency
% calculation
if strcmpi(s.costfunction,'torque') || strcmpi(s.costfunction,'covariant acceleration') || strcmpi(s.costfunction,'acceleration coord') || strcmpi(s.costfunction,'power quality')
    % With cost as time period, period of the gait is the cost to execute
    % the gait at unit torque squared to the 1/4th power
    totalstroke = cost^(1/4);
else
    % Drag systems have cost as path length, so no modification is needed
    totalstroke = cost;
end

% If efficiency is negative, reversing the order of points so that
% efficiency is positive
if lineint<0
    lineint= -lineint;
    y=flip(y);
    invert=1;
else
    invert=0;
end

%% Preliminaries for gradient calculation
ccf = ccf_interpolator(y,s,n,dimension,direction);
[metric,metricgrad] = metric_interpolator(y,s,n,dimension);

%% Jacobianstroke is the gradient of cost. 
%Contrigrad is the contribution to the gradient ue to the metric changing

jacobian = struct();

switch s.costfunction
    case {'pathlength metric','pathlength coord','pathlength metric2'}
        % Get the gradient of cost based on drag-dominated system
        jacobian.stroke = jacobianstrokecalculator(y,n,dimension,metric,metricgrad);
    case {'torque','covariant acceleration','acceleration coord','power quality'}
        % Get the gradient of cost based on inertia-dominated system
        inertia_cost_grad = inertia_cost_gradient(s,n,coeff,T,p,'discrete');
        
        % With cost as time period to execute the gait, the gradient of
        % cost for inertial systems is the gradient of cost with period 1
        % divided by (4*T^3)
        inertia_cost_grad = inertia_cost_grad./(4*totalstroke^3);
    otherwise
        error('Unexpected system type at cost gradient calculation!')
end

%% Jacobiandisp is the gradient of displacement.
% jacobiandispcalculator3 is the function that calculates the gradient of 
% displacement for the ith point. It's input arguments are the coordinates of 
% the (i-1)th, ith and (i+1)th point, CCF value at point i and the dimension of     
% the shape space (dimension)


jacobian.disp = zeros(n,dimension);
for i=2:1:n-1
    jacobian.disp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),ccf(i,:),dimension);
end
jacobian.disp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),ccf(1,:),dimension);
jacobian.disp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),ccf(n,:),dimension);

%% Jacobianeqi is the concentration gradient. 
%It is the term that keeps points eqi distant from each other and prevents crossover of gait.

jacobian.eqi = jacobianeqicalculator(y,n,dimension,metric);

%% Jacobina.repuls is the repulsive vector field for a joint limit.
jacobian.repuls = zeros(n,dimension);
if ~isempty(lb)
    lb = reshape(lb,[n+1 dimension]);
    lb = lb(1:n,:);
    jacobian.repuls = jacobian.repuls + 5*sech(10*(y-lb)).^2;
end

if ~isempty(ub)
    ub = reshape(ub,[n+1 dimension]);
    ub = ub(1:n,:);
    jacobian.repuls = jacobian.repuls + 5*sech(10*(ub-y)).^2;
end

%% properly ordering gradients depending on whether lineint was negative or positive
if invert
    y=flip(y);
    jacobian.disp=flip(jacobian.disp);
    jacobian.eqi=flip(jacobian.eqi);
    jacobian.stroke=flip(jacobian.stroke);
end

end

