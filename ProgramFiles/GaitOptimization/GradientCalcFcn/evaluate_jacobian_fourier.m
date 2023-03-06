function [jacobfourier,net_disp_opt,totalstroke] = evaluate_jacobian_fourier(y,s,n,dimension,direction,lb,ub,~)
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
% jacobfourier: the struct which contains the gradient of functions.
%   (jacobfourier.disp, jacobfourier.stroke,
%   jacobfourier.eqi, jacobfourier.eff)
% lineint : the displacement for the gait
% totalstroke : the cost for the gait
%%%%%%%%%%%%%


%% Check if there is the inequality condition (e.g. joint limit)
if nargin < 6
    lb = [];
    
    ub = [];
end

%% Obtaining points from fourier coefficients
% The first step is to obtain a direct transcription of the gait from the
% fourier series parametrization. y1 is the matrix of coordinate values of
% the points forming the gait.
coeff=y;
y1 = path_from_fourier(y,n,dimension);
y1 = y1(1:end-1,:); % Remove the last points because path_from_fourier returns self-connected gait

%% Calculating cost and displacement per gait
w1 = y(end,1); % Frequency of Fourier transform
w2 = y(end,2);
% Assign a time period for executing the gait
T = 2*pi/w1;

% Define phi_def = [alpha1, alpha2] as a function of time t such that the
% array returns the shape variables given by the fourier coefficients at a
% time t
p = makeGait(y);
        
% % Uncomment this section to verify that the shape variables and derivatives
% % have been derived appropriately
% valid = verify_shape_equations(p);
% valid_M = verify_mass_derivative(s);
% validate_cost_gradient(s,n,y,T,p);

% Reassignment of point transcription to variable y for compatibility with
% old version
clear y
y=y1;
clear y1;

% Calculate displacement, cost and efficiency of a gait
% Note that, for inertial cost, cost is returned as the integral of torque
% squared, while for drag-based systems, cost is the path length of the
% gait
[~, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');

% The coordinate of gradient is on the exponential coordinate, 
% but the net displacement derived by line integral is on the original
% coordinate. Take the matrix logarithms on the net displacement.
net_disp_opt = se2log(net_disp_opt);

% displacement produced in the chosen direction produced on executing 
% the gait measured in the optimal coordinates 
lineint=net_disp_opt(direction);

%% Checking if the system is inertia-dominated or drag-dominated
switch s.costfunction
    case {'pathlength metric','pathlength coord','pathlength metric2'}
        % drag-dominated system
        systemtype = 0;
    case {'torque','covariant acceleration','acceleration coord','power quality'}
        % inertia-dominated system
        systemtype = 1;
    otherwise
        error('Unexpected system type at cost gradient calculation!')
end

% Assign value for totalstroke, i.e. the cost metric used for efficiency
% calculation
if systemtype
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

%% changey/dcoeff tells us how much each point moves when a fourier series variable is changed
chy=chy_generator(coeff,n,dimension);


%% Jacobian.stroke is the gradient of cost. 
%Contrigrad is the contribution to the gradient ue to the metric changing

jacobian = struct();
if systemtype
    % Get the gradient of cost based on inertia-dominated system
    inertia_cost_grad = inertia_cost_gradient(s,n,coeff,T,p,'discrete');

    % With cost as time period to execute the gait, the gradient of
    % cost for inertial systems is the gradient of cost with period 1
    % divided by (4*T^3)
    inertia_cost_grad = inertia_cost_grad./(4*totalstroke^3);
else
    % Get the gradient of cost based on drag-dominated system
    jacobian.stroke = jacobianstrokecalculator(y,n,dimension,metric,metricgrad);
end

%% Jacobian.disp is the gradient of displacement.
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

%% Jacobian.eqi is the concentration gradient. 
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


%% properly ordering gradients depending on wether lineint was negative or positive
if invert
    jacobian.disp=flip(jacobian.disp);
    jacobian.eqi=flip(jacobian.eqi);
    jacobian.stroke=flip(jacobian.stroke);
    jacobian.repuls=flip(jacobian.repuls);
end


%% fourier series version of all gradients

% We then obtain gradients in a fourier series parametrization by
% projecting the gradients from the direct transcription space onto the
% fourier coefficient space
jacobfourier=struct();
jacobfourier.disp = zeros(length(coeff)-1,dimension);
jacobfourier.stroke = zeros(length(coeff)-1,dimension);
jacobfourier.eqi = zeros(length(coeff)-1,dimension);
jacobfourier.repuls = zeros(length(coeff)-1,dimension);

for i=1:1:dimension
    for j=1:1:length(coeff)-1
        jacobfourier.disp(j,i)=chy{i}(j,:)*jacobian.disp(:,i);
        if systemtype == 0
            jacobfourier.stroke(j,i)=chy{i}(j,:)*jacobian.stroke(:,i);
        end
        jacobfourier.eqi(j,i)=chy{i}(j,:)*jacobian.eqi(:,i);
        jacobfourier.repuls(j,i)=chy{i}(j,:)*jacobian.repuls(:,i);
    end
end

if systemtype == 1

    % Inertia cost gradient is already in terms of the fourier coefficients
    jacobfourier.stroke = inertia_cost_grad;

    % Add zero terms to the gradient of displacement and jacobianeqi for
    % frequency terms, since the inertia gradient includes those terms
    jacobfourier.disp = [jacobfourier.disp;zeros(1,dimension)];
    jacobfourier.eqi = [jacobfourier.eqi;zeros(1,dimension)];    % Calculate the total gradient of efficiency

    % jacobianeqifourier commented out at Hossein's suggestion - don't
    % necessarily want the points to all be equally spaced
    jacobfourier.eff = jacobfourier.disp/totalstroke-lineint*jacobfourier.stroke/totalstroke^2;%+jacobianeqifourier;

    % reset the gradient of the frequency terms to be zero so they aren't
    % changed
    jacobfourier.eff(end,:) = zeros(1,dimension);
else    
    jacobfourier.eff = jacobfourier.disp-...
        lineint*jacobfourier.stroke/totalstroke+jacobfourier.eqi;
end

end

