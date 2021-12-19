function [f,g]=solvedifffmincon(y,s,n,dimension,direction,~,~,~)%,lb,ub,writerObj)
%%%%%%%%%%%%%
% This function calculates efficiency (or displacement, if
% that is the objective function) and its gradient with respect to the coefficients obtained
% by the fourier series parametrization
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
% lb: Lower bound of shape variables for each point which is obtained from the grid inside 
%   which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside 
%   which an optimal gait is desired
% 
% Outputs:
% 
% f: Objective function value (This is negative of efficiency by default, can be
%   changed to displacement)
% g: Gradient of the objective function
%%%%%%%%%%%%%

%% Obtaining points from fourier coefficients
% The first step is to obtain a direct transcription of the gait from the
% fourier series parametrization. y1 is the matrix of coordinate values of
% the points forming the gait.
afactor=0.001;
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

global bestCost bestDisp bestEff;

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
    ytemp=y;
    for i=1:n
        y(i,:)=ytemp(n+1-i,:);
    end
    invert=1;
else
    invert=0;
end

%% Preliminaries for gradient calculation
% Preallocating memory for variables which we will need in further
% calculation 
yvalues=cell(n,dimension); % Cell representation of the coordinates of all points forming the gait
interpstateccf=cell(1,dimension); % Variable which will store the ccf function grid used for interpolation
interpmetricgrid=cell(1,dimension);  % Variable which will store the metric grid used for interpolation
ccf=zeros(n,dimension*(dimension-1)/2); % Variable which will store ccf function at each point
metric1=zeros(n,dimension,dimension); % Variable which will store metric at each point in the form of a matrix
metric = repmat({zeros(dimension)},[n 1]); % Variable which stores the metric at each point in the form of a 2x2 matrix
metricgrad1=zeros(n,dimension,dimension,dimension); % Variable which will store gradient of metric at each point in the form of a matrix
metricgrad = repmat({zeros(dimension)},[n,dimension]); % Variable which will store gradient of metric at each point in the form of a matrix

% Interpolation to calculate all the variables needed for gradient
% calculation
for i=1:1:n
    for j=1:1:dimension
        yvalues{i,j}=y(i,j);
    end
end

y_for_interp = mat2cell(y,size(y,1),ones(1,size(y,2)));

for j=1:1:dimension
    interpstateccf{j}=s.grid.eval{j,1};
    interpmetricgrid{j}=s.grid.metric_eval{j,1};
end

for j=1:dimension*(dimension-1)/2
    ccf(:,j)=interpn(interpstateccf{:},s.DA_optimized{direction,j},y_for_interp{:},'spline');
end

for j=1:1:dimension
    for k=1:1:dimension
        metric1(:,j,k)=interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y_for_interp{:},'spline');
    end
end

if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'acceleration coord')
    for i=1:n
       metric{i}=eye(dimension);
    end
elseif strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
    for i=1:n
        for j=1:1:dimension
           for k=1:1:dimension
               metric{i}(j,k)=metric1(i,j,k);
           end
        end
    end
end
if strcmpi(s.costfunction,'pathlength metric2')
    for i = 1:n
        metric{i} = metric{i}*metric{i};
    end
end

if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'acceleration coord')
    for i=1:dimension
        for j=1:n
            metricgrad{j,i} = zeros(dimension);
        end
    end
elseif strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
    y2 = zeros(size(y));
    y1 = y2;
    for l=1:1:dimension
        for m=1:1:dimension
            if m==l
               y2(:,m)=y(:,m)+afactor*ones(length(y),1);
               y1(:,m)=y(:,m)-afactor*ones(length(y),1);
            else
               y2(:,m)=y(:,m);
               y1(:,m)=y(:,m);
            end
        end
        y2_for_interp = mat2cell(y2,size(y,1),ones(1,size(y,2)));
        y1_for_interp = mat2cell(y1,size(y,1),ones(1,size(y,2)));
        for j=1:1:dimension
            for k=1:1:dimension
                metricgrad1(:,l,j,k)=(interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y2_for_interp{:},'spline')...
                    -interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y1_for_interp{:},'spline'))/(2*afactor);
            end
        end
        for i=1:n
            for j=1:1:dimension
                for k=1:1:dimension
                    metricgrad{i,l}(j,k)=metricgrad1(i,l,j,k);
                end
            end
        end
    end
end
if strcmpi(s.costfunction,'pathlength metric2')
    for l = 1:dimension
        for i = 1:n
            metricgrad{i,l} = metricgrad{i,l}*metric{i}+metric{i}*metricgrad{i,l};
        end
    end
end

%% changey/dcoeff tells us how much each point moves when a fourier series variable is changed
% chy is a cell with as many entries as the dimension of the shape space
% ith element of chy is a matrix where the (j,k)th entry tells us the change in the ith coordinate
% of the kth point of the gait resulting from a unit change in the jth
% fourier coefficient corresponding to the ith dimension of the shape space

chy=cell(dimension,1);
% Create vector of time values at which to evaluate points of gait
t = linspace(0,T,n);
for i=1:1:dimension
    for j=1:1:n
        chy{i}(:,j)=[1;cos(t(j)*coeff(end,i));sin(t(j)*coeff(end,i));cos(2*t(j)*coeff(end,i));sin(2*t(j)*coeff(end,i));cos(3*t(j)*coeff(end,i));sin(3*t(j)*coeff(end,i));cos(4*t(j)*coeff(end,i));sin(4*t(j)*coeff(end,i))];%cos(5*t(j)*coeff(end,i));sin(5*t(j)*coeff(end,i))];%;cos(6*t(j)*coeff(end,i));sin(6*t(j)*coeff(end,i))];%
    end
end

%% Jacobianstroke is the gradient of cost. 
%Contrigrad is the contribution to the gradient due to the metric changing
switch s.costfunction
    case {'pathlength metric','pathlength coord','pathlength metric2'}
        % Get the gradient of cost based on drag-dominated system
        jacobianstroke = jacobianstrokecalculator(y,n,dimension,metric,metricgrad);
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

jacobiandisp = zeros(n,dimension);
for i=2:1:n-1
    jacobiandisp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),ccf(i,:),dimension);
end
jacobiandisp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),ccf(1,:),dimension);
jacobiandisp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),ccf(n,:),dimension);

%% Jacobianeqi is the concentration gradient. 
%It is the term that keeps points eqi distant from each other and prevents crossover of gait.

jacobianeqi = jacobianeqicalculator(y,n,dimension,metric);

%% properly ordering gradients depending on wether lineint was negative or positive
if invert
        jacobiandisptemp=jacobiandisp;
        if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
                jacobianstroketemp=jacobianstroke;
        end
        jacobianeqitemp=jacobianeqi;
    for i=1:1:n
        jacobiandisp(i,:)=jacobiandisptemp(n+1-i,:);
        if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
            jacobianstroke(i,:)=jacobianstroketemp(n+1-i,:);
        end
        jacobianeqi(i,:)=jacobianeqitemp(n+1-i,:);
    
    end
end

%% fourier series version of all gradients

% We then obtain gradients in a fourier series parametrization by
% projecting the gradients from the direct transcription space onto the
% fourier coefficient space
for i=1:1:dimension
    for j=1:1:9 
        jacobiandispfourier(j,i)=chy{i}(j,:)*jacobiandisp(:,i);
        if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
            jacobianstrokefourier(j,i)=chy{i}(j,:)*jacobianstroke(:,i);
        end
        jacobianeqifourier(j,i)=chy{i}(j,:)*jacobianeqi(:,i);
    end
end
if strcmpi(s.costfunction,'torque') || strcmpi(s.costfunction,'covariant acceleration') || strcmpi(s.costfunction,'acceleration coord') || strcmpi(s.costfunction,'power quality')
    % Inertia cost gradient is already in terms of the fourier coefficients
    jacobianstrokefourier = inertia_cost_grad;
    % Add zero terms to the gradient of displacement and jacobianeqi for
    % frequency terms, since the inertia gradient includes those terms
    jacobiandispfourier = [jacobiandispfourier;zeros(1,dimension)];
    jacobianeqifourier = [jacobianeqifourier;zeros(1,dimension)];
    % Calculate the total gradient of efficiency
    % jacobianeqifourier commented out at Hossein's suggestion - don't
    % necessarily want the points to all be equally spaced
    totaljacobianfourier = jacobiandispfourier/totalstroke-lineint*jacobianstrokefourier/totalstroke^2;%+jacobianeqifourier;
    % reset the gradient of the frequency terms to be zero so they aren't
    % changed
    totaljacobianfourier(end,:) = zeros(1,dimension);
else
    totaljacobianfourier = jacobiandispfourier/totalstroke-lineint*jacobianstrokefourier/totalstroke^2+jacobianeqifourier;
end

%% minimizing negative of efficiency(or displacement)
f=-lineint/(totalstroke); % Optimizing for displacement over cost
%     f = -lineint; % Optimizing for displacement only

if abs(f) > bestEff
    bestEff = abs(f);
    bestDisp = abs(lineint);
    bestCost = abs(totalstroke);
end

if nargout>1
    if strcmpi(s.costfunction,'torque') || strcmpi(s.costfunction,'covariant acceleration') || strcmpi(s.costfunction,'acceleration coord') || strcmpi(s.costfunction,'power quality')
        % Return the gradient of efficiency as previously calculated for
        % inertia systems
        g = -totaljacobianfourier; % Optimizing for displacement over cost
%         g = -jacobiandispfourier; % Optimizing for displacement only
    else
        % Return the gradient of efficiency plus row of zeros for frequency
        % terms for drag systems
        g=[-totaljacobianfourier;zeros(1,dimension)]; % Optimizing for displacement over cost
%         g = [-jacobiandispfourier;zeros(1,dimension)]; % Optimizing for displacement only
    end
end

end