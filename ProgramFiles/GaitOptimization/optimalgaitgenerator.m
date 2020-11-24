function y=optimalgaitgenerator(s,dimension,npoints,a,lb,ub,stretch,direction,costfunction,handles)
%%%%%%%%%%%%%%
% This function takes an input gait and runs fmincon to find the neareast locally 
% optimal gait

%Inputs:
%
%s: System file which contains the connection vector field, CCF's and
%   metric data
%dimension: Indicates the number of shape variables of the system
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
% a: n x m array, for n waypoints in m dimensions
% lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% direction: Direction to optimize travel: 1=x,2=y,3=theta
% costfunction: costfunction type to optimize for
%           Options: Path length (metric), torque squared (metric),
%           covariant acceleration (metric), path length (coordinate),
%           acceleration (coordinate)
% 
% 
% Outputs: 
% 
% y: Matrix whose values indicate coordinates of points which form the optimal gait
%%%%%%%%%%%%

% n=npoints;
% for i=1:1:npoints
%     a1(1,i)=1*cos(2*pi*(i-1)/n);
%     a2(1,i)=1*cos(2*pi*(i-1)/n+pi/2);
% end

% n=npoints;
% P1(:,1)=a1(1,1:n)';
% P1(:,2)=a2(1,1:n)';
% P1(end+1,:) = P1(1,:); % Close the loop on the gait

% For minimal refactoring, mapping a -> P1
P1 = a;

% % Close the loop of the gait if necessary
% if P1(end,:) ~= P1(1,:)
%     P1(end+1,:) = P1(1,:);
% end


%% Finding fourier coeffecients.
% The first step is to go from a direct transcription of the initial gait
% to a fourier based parametrization. 
% fa is a cell where the ith element contains the coefficients for the fourier parametrization along the ith direction 

% Time period of gait is 1 second in order to handle calculations performed
% for inertial gait optimization
t = linspace(0,1,size(P1,1));

fa=cell(dimension);
% The bounds ensure the fourier series terms have the right period
 % Bound only the frequency to be 2*pi, such that period = 1
options = fitoptions('fourier4');
options.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 2*pi];
options.Upper = -[-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -2*pi];

for i=1:1:dimension
    fa{i}=fit(t',P1(:,i),'fourier4',options);
end

%% The next step is to setup the fmincon call. 
% y0 is the marix of all fourier series coefficients that describe the
%   initial gait
% nonlcon is the function that imposes the constraints that all the points
%   stay inside the grid
% outfun is the function that updates the gait on the GUI after every iteration 

A=[];
b=[];
Aeq=[];
beq=[];

nu={'a0';'a1';'b1';'a2';'b2';'a3';'b3';'a4';'b4';'w'};%
 lb1=[];
 ub1=[];

 y0 = zeros(length(nu),dimension);
for i=1:dimension
    for j=1:length(nu)
        y0(j,i)=fa{i}.(nu{j});
    end
end

writerObj = [];
% % Uncomment this section if you'd like to record a video of the
% % optimizer's steps in the shape space
% writerObj = VideoWriter('cost_as_time_period_circle_start.mp4','MPEG-4');
% writerObj.FrameRate = 5;
% % set the seconds per image
% % open the video writer
% open(writerObj);
% figure(5);
% subplot(1,2,1)
% contour(s.grid.eval{1},s.grid.eval{2},s.DA_optimized{1},'LineWidth',1.5)
% axis square
% hold on

s.costfunction = costfunction;
if strcmpi(costfunction,'torque')
    s.metric = @(alpha1,alpha2)eye(2);
end

%Suppress warning for annoying thing in jacobianeqicalculator
warning('off','MATLAB:sqrtm:SingularMatrix');

try
 options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','iter','Algorithm','sqp','CheckGradients',false,'FiniteDifferenceType','central','MaxIter',4000,'MaxFunEvals',20000,'TolCon',10^-2,'OutputFcn', @(y,optimValues,state) outfun(y,optimValues,state,stretch,s,handles));
catch
   error('This code requires the global optimization toolbox to run') 
end
 [yf, ~,~,~]=fmincon(@(y) solvedifffmincon(y,s,npoints,dimension,direction,lb,ub,writerObj),y0,A,b,Aeq,beq,lb1,ub1,@(y) nonlcon(y,s,npoints,dimension,lb,ub),options);

% % Uncomment this if you uncommented the section above so that the video
% % writer object is closed appropriately.
% close(writerObj);

%% Getting point position values from the result of fmincon
% This section helps us go back to a direct transcription parametrization
% of the optimal gait from a fourier series parametrization. y is a column vector
% that contains coordinates of all points forming the optimized gait

y1 = path_from_fourier(yf,npoints,dimension);
% path_from_fourier returns a self-connected gait, so remove the last point
% to give what optimalgaitgenerator expects to return
y1 = y1(1:end-1,:);
y=y1(:);

%% Uncomment for plotting the optimized gait. Potentially useful while debugging.
% for i=1:n
%     xf(i)=y(i);
%     yf(i)=y(n+i);
% end
% % 
% % for i=1:n`
% %     xf(i)=P1(i,1);
% %     yf(i)=P1(i,2);
% % end
% % 
% figure(11)
% hold on
% plot(xf,yf,'-o')

end

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
p.phi_def = @(t) [y(1,1)+y(2,1)*cos(w1*t)+y(3,1)*sin(w1*t)+y(4,1)*cos(2*w1*t)+...
            +y(5,1)*sin(2*w1*t)+y(6,1)*cos(3*w1*t)+y(7,1)*sin(3*w1*t)+...
            +y(8,1)*cos(4*w1*t)+y(9,1)*sin(4*w1*t),...
            y(1,2)+y(2,2)*cos(w2*t)+y(3,2)*sin(w2*t)+y(4,2)*cos(2*w2*t)+...
            +y(5,2)*sin(2*w2*t)+y(6,2)*cos(3*w2*t)+y(7,2)*sin(3*w2*t)+...
            +y(8,2)*cos(4*w2*t)+y(9,2)*sin(4*w2*t)];

% Similarly for the shape velocities
p.dphi_def = @(t) [-(w1)*y(2,1)*sin(w1*t)+(w1)*y(3,1)*cos(w1*t)-2*(w1)*y(4,1)*sin(2*w1*t)+...
            +2*(w1)*y(5,1)*cos(2*w1*t)-3*(w1)*y(6,1)*sin(3*w1*t)+3*(w1)*y(7,1)*cos(3*w1*t)+...
            -4*(w1)*y(8,1)*sin(4*w1*t)+4*(w1)*y(9,1)*cos(4*w1*t),...
            -(w2)*y(2,2)*sin(w2*t)+(w2)*y(3,2)*cos(w2*t)-2*(w2)*y(4,2)*sin(2*w2*t)+...
            +2*(w2)*y(5,2)*cos(2*w2*t)-3*(w2)*y(6,2)*sin(3*w2*t)+3*(w2)*y(7,2)*cos(3*w2*t)+...
            -4*(w2)*y(8,2)*sin(4*w2*t)+4*(w2)*y(9,2)*cos(4*w2*t)];

% Similarly for the shape accelerations
p.ddphi_def = @(t) [-(w1)^2*y(2,1)*cos(w1*t)-(w1)^2*y(3,1)*sin(w1*t)-4*(w1)^2*y(4,1)*cos(2*w1*t)+...
            -4*(w1)^2*y(5,1)*sin(2*w1*t)-9*(w1)^2*y(6,1)*cos(3*w1*t)-9*(w1)^2*y(7,1)*sin(3*w1*t)+...
            -16*(w1)^2*y(8,1)*cos(4*w1*t)-16*(w1)^2*y(9,1)*sin(4*w1*t),...
            -(w2)^2*y(2,2)*cos(w2*t)-(w2)^2*y(3,2)*sin(w2*t)-4*(w2)^2*y(4,2)*cos(2*w2*t)+...
            -4*(w2)^2*y(5,2)*sin(2*w2*t)-9*(w2)^2*y(6,2)*cos(3*w2*t)-9*(w2)^2*y(7,2)*sin(3*w2*t)+...
            -16*(w2)^2*y(8,2)*cos(4*w2*t)-16*(w2)^2*y(9,2)*sin(4*w2*t)];
        
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
[~, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','ODE');
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

for j=1:1:dimension
    interpstateccf{j}=s.grid.eval{j,1};
    interpmetricgrid{j}=s.grid.metric_eval{j,1};
end

for j=1:dimension*(dimension-1)/2
    ccf(:,j)=interpn(interpstateccf{:},s.DA_optimized{direction,j},y(:,1),y(:,2),'spline');
end

for j=1:1:dimension
    for k=1:1:dimension
        metric1(:,j,k)=interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y(:,1),y(:,2),'spline');
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
        for j=1:1:dimension
            for k=1:1:dimension
                metricgrad1(:,l,j,k)=(interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y2(:,1),y2(:,2),'spline')-interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y1(:,1),y1(:,2),'spline'))/(2*afactor);
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


%% Variables useful while debugging for flaws in gradient calculations

% for i=1:n
%     for j=1:1:dimension
%         totaljacobianc(i,j)=chy{j}(:,i)'*totaljacobianfourier(:,j);
%         jacobiandispc(i,j)=chy{j}(:,i)'*jacobiandispfourier(:,j);
%         jacobianstrokec(i,j)=chy{j}(:,i)'*jacobianstrokefourier(:,j);
%     end
% end
%         totaljacobianctemp=totaljacobianc;
%         jacobiandispctemp=jacobiandispc;
%         jacobianstrokectemp=jacobianstrokec;
% for i=1:1:n
%     jacobiandispc(i,:)=jacobiandispctemp(n+1-i,:);
%     jacobianstrokec(i,:)=jacobianstrokectemp(n+1-i,:);
%     totaljacobianc(i,:)=totaljacobianctemp(n+1-i,:);
% end

%% minimizing negative of efficiency(or displacement)
 f=-lineint/(totalstroke); % Optimizing for displacement over cost
%     f = -lineint; % Optimizing for displacement only
%  % Uncomment these lines if you'd like to print out the displacement
%  % (lineint) and cost (totalstroke) while the program is optimizing
%  lineint
%  totalstroke
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

%% Debugging and plotting
% This section was written up for the sole purpose of helping with
% debugging flaws in the optimizer. Appropriate sections can be uncommented
% to plot how the different gradients look during the optimization process.

% for i=1:n
%     G(i)=y(i,1);
%     H(i)=y(i,2);
% %     P(i)=y(i,3);
% %     I(i)=1*jacobianstroketemp(i,1);
% %     J(i)=1*jacobianstroketemp(i,2);
% %     I(i)=1*jacobianstroke(i,1);
% %     J(i)=1*jacobianstroke(i,2);
% %     N(i)=1*jacobianstroke(i,3);
% %     K(i)=jacobiandisptemp(i,1);
%     K(i)=jacobiandisp(i,1);
% %     K1(i) = jacobiandispc(i,1);
% %     K1(i)=jacobianforward(i,1);
% %     L(i)=jacobiandisptemp(i,2);
%     L(i)=jacobiandisp(i,2);
%     L1(i)=jacobiandispc(i,2);
%     L1(i)=jacobianforward(i,2);
%     Q(i)=jacobiandisp(i,3);
%     B(i)=totaljacobianc(i,1);
%     B1(i)=totaljacobian(i,1);
%     C(i)=totaljacobianc(i,2);
%     C1(i)=totaljacobian(i,2);
%     D(i)=totaljacobian(i,3);
%     O(i)=jacobianeqi(i,1);
%     S(i)=jacobianeqi(i,2);
%     R(i)=jacobianeqi(i,3);
%     ccfx(i)=ccf(i,1);
%     ccfy(i)=ccf(i,2); 
%     ccfz(i)=ccf(i,3);
%     ccfcrosscheck(i)=ccf(i,:)*jacobiandisp(i,:)';
% end
% % % % 
% % % % 
% clf(figure(6)) %%jacobiandisp
% figure(6)
% scale=0;
% scale1=1;
% quiver(G,H,K,L,scale1)
% hold on
% % quiver(G,H,K1,L1,scale1)
% % quiver(G,H,ccfx,ccfy,scale1)
% plot(G,H)
% axis equal
% hold off
% 
% % clf(figure(9)) %%jacobianforward
% % figure(9)
% % scale=0;
% % scale1=1;
% % quiver(G,H,K1,L1,scale1)
% % hold on
% % % quiver(G,H,P,ccfx,ccfy,ccfz,scale1)
% % plot(G,H)
% % axis equal
% % hold off
% % 
% 
% 
% % % 
% clf(figure(7)) %%jacobianstroke
% figure(7)
% scale=0;
% quiver(G,H,K,L,scale)
% quiver(G,H,I,J,scale1)
% hold on
% plot(G,H)
% axis equal
% hold off
% % % % 
% clf(figure(8)) %%% totaljacobian
% figure(8)
% scale=0;
%  quiver(G,H,B,C,scale1)
% hold on
% quiver(G,H,B1,C1,scale1)
% 
% plot(G,H)
% axis equal
% hold off
% % 
% % clf(figure(3)) %%%jacobianeqi
% % figure(3)
% % scale1=1;
% % quiver3(G,H,P,O,S,R,scale1)
% % hold on
% % plot3(G,H,P)
% % axis equal
% % hold off
% % 
% % 
% % % figure(5)
% % % scale=0;
% % % quiver(G,H,B1,C1,scale1)
% % % hold on
% % % plot(G,H)
% % % axis equal
% % % hold off
% % min(abs(totaljacobian))
% % max(abs(totaljacobian))
% % pause(0.1)
% 
% % 
% % for i=1:1:length(ccf)
% %     ccfcheckper(i)=(ccfcheck(i)-ccf(i))/ccf(i);
% % end
% 

% % Uncomment this section if you'd like to plot the optimizer's progress
% % more frequently (i.e. not just the steps that worked) and save a video of
% % the frames
% figure(5)
% subplot(1,2,1)
% plot(y(:,1),y(:,2),'r','LineWidth',2)
% title({['Disp: ',num2str(lineint)],['Cost: ',num2str(totalstroke)],['Period from freq: ',num2str(2*pi/w1)]})
% 
%  pause(0.05)
%  frame = getframe(gcf);
%  writeVideo(writerObj,frame);
end

function jacobianstroke = jacobianstrokecalculator(y,n,dimension,metric,metricgrad)
    % Calculates the gradient of cost for drag dominated systems.
    % Inputs:
    %   y: matrix containing points that transcribe the gait
    %   n: number of points in gait transcription
    %   dimension: number of shape variables
    %   metric: Riemannian metric
    %   metricgrad: Gradient of Riemannian metric
    
    %l is the vector containing metric weighted distances between neighbouring
    %points
    l = zeros(1,n);
    for i=1:(numel(l)-1)
        l(i)=sqrt((y(i+1,:)-y(i,:))*((metric{i}+metric{i+1})/2)*(y(i+1,:)-y(i,:))');
    end
    l(end)=sqrt((y(1,:)-y(n,:))*((metric{n}+metric{1})/2)*(y(1,:)-y(n,:))');
    
    delp = cell(1,n);
    for i=1:(numel(delp)-1)
        delp{i}=y(i+1,:)-y(i,:); % delp{i} is the vector joining the (i+1)th point to the ith point 
    end
    delp{end}=y(1,:)-y(n,:);

    jacobianstroke = zeros(n,dimension);
    contrigrad=zeros(n,dimension);
    for i=2:n-1
        for j=1:dimension
            %Contrigrad is the contribution to the gradient due to the metric changing
            contrigrad(i,j)=0.5*delp{i}*metricgrad{i,j}*delp{i}'/(2*l(i))+0.5*delp{i-1}*metricgrad{i,j}*delp{i-1}'/(2*l(i-1)); 
        end
        % Total gradient is the result of distance changing due to movement of point and the metric changing due to movement of the point
        jacobianstroke(i,:)=(-(((metric{i}+metric{i+1})/2)*delp{i}')'-(delp{i}*((metric{i}+metric{i+1})/2)))/(2*l(i))+...
            +((((metric{i-1}+metric{i})/2)*delp{i-1}')'+(delp{i-1}*((metric{i}+metric{i-1})/2)))/(2*l(i-1))+contrigrad(i,:); 
    end

    % Calculation for the 1st point and last point have to be done outside the
    % loop as the (i+1)th point for the last point is the first point and
    % (i-1)th point for the first point is the last point
    for j=1:dimension
        contrigrad(1,j)=0.5*delp{1}*metricgrad{1,j}*delp{1}'/(2*l(1))+0.5*delp{n}*metricgrad{1,j}*delp{n}'/(2*l(n));
    end
    jacobianstroke(1,:)=(-(((metric{1}+metric{2})/2)*delp{1}')'-(delp{1}*((metric{1}+metric{2})/2)))/(2*l(1))+...
        +((((metric{n}+metric{1})/2)*delp{n}')'+(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+contrigrad(1,:);

    for j=1:dimension
        contrigrad(n,j)=0.5*delp{n}*metricgrad{n,j}*delp{n}'/(2*l(n))+0.5*delp{n-1}*metricgrad{n,j}*delp{n-1}'/(2*l(n-1));
    end
    jacobianstroke(n,:)=(-(((metric{n}+metric{1})/2)*delp{n}')'-(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+...
        +((((metric{n}+metric{n-1})/2)*delp{n-1}')'+(delp{n-1}*((metric{n}+metric{n-1})/2)))/(2*l(n-1))+contrigrad(n,:);
end

function jacobianeqi = jacobianeqicalculator(y,n,dimension,metric)
    % Calculates the gradient of the force driving the points along the
    % gait to be equally spaced.
    % Inputs:
    %   y: matrix containing points that transcribe the gait
    %   n: number of points in gait transcription
    %   dimension: number of shape variables
    %   metric: Riemannian metric
    
    jacobianeqi = zeros(n,dimension);
    
    %l is the vector containing metric weighted distances between neighbouring
    %points
    l = zeros(1,n);
    for i=1:(numel(l)-1)
        l(i)=sqrt((y(i+1,:)-y(i,:))*((metric{i}+metric{i+1})/2)*(y(i+1,:)-y(i,:))');
    end
    l(end)=sqrt((y(1,:)-y(n,:))*((metric{n}+metric{1})/2)*(y(1,:)-y(n,:))');
    
    for i=2:n-1
        len=sqrt((y(i+1,:)-y(i-1,:))*((metric{i-1}+metric{i+1})/2)*(y(i+1,:)-y(i-1,:))'); % metric weighted length between point (i-1) and (i+1)
        midpoint=y(i-1,:)+((y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2))/2; % location of midpoint of the line segment joining point (i-1) and (i+1)
        betacos=(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*((y(i,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i})/2))'/(l(i-1)*len);
        xhat=y(i-1,:)+(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*l(i-1)*betacos/len; %projection of ith point onto the line joining the (i-1)th and (i+1)th points
        jacobianeqi(i,:)=midpoint-xhat; % gradient of the ith point is equal to the difference between the midpoint and the projection of ith point
    end

    len=sqrt((y(2,:)-y(n,:))*((metric{2}+metric{n})/2)*(y(2,:)-y(n,:))');
    midpoint=y(n,:)+((y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2))/2;
    betacos=(y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2)*((y(1,:)-y(n,:))*sqrtm((metric{n}+metric{1})/2))'/(l(n)*len);
    xhat=y(n,:)+(y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2)*l(n)*betacos/len;
    jacobianeqi(1,:)=midpoint-xhat;

    len=sqrt((y(1,:)-y(n-1,:))*((metric{1}+metric{n-1})/2)*(y(1,:)-y(n-1,:))');
    midpoint=y(n-1,:)+((y(1,:)-y(n-1,:))*sqrtm((metric{1}+metric{n-1})/2))/2;
    betacos=(y(1,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{1})/2)*((y(n,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{n})/2))'/(l(n-1)*len);
    xhat=y(n-1,:)+(y(1,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{1})/2)*l(n-1)*betacos/len;
    jacobianeqi(n,:)=midpoint-xhat;
end

function a=jacobiandispcalculator3(p1,p2,p3,ccf,dimension)
%%%%%%%%%
%
% jacobiandispcalculator3 is the function that calculates the gradient of 
% displacement for the ith point. 
% Its input arguments are the coordinates of the (i-1)th, ith and (i+1)th point,     
% CCF value at point i(ccf) and the dimension of the shape space (dimension)
%
%%%%%%%%%

l1=0; % variable for calculating length of the line segment joining the (i-1)th point with the (i+1)th point
base = zeros(1,dimension);
for i=1:numel(base)
    l1=l1+(p1(i)-p3(i))^2;
    base(1,i)=p3(i)-p1(i); % vector connecting the (i-1)th point and (i+1)th point  
end
%l=sqrt(l1); % length of the line segment joining the (i-1)th point with the (i+1)th point

jacobian = zeros(1,dimension);
for i=1:dimension
%    jacobian(1,i)=0;
    perp1=zeros(1,dimension);
    perp1(i)=1;
    %parcomp=base*perp1'/norm(base);
    %perp1-parcomp*base/norm(base);  %%recheck again
    perp=perp1;% Unit vector along the ith direction
    % The for loop below calculates the gradient along the ith direction by
    % treating the CCF as 2 forms. A specific (j,k) represents a component of the 2 form 
    for j=1:dimension-1 
        for k=1:dimension-j
            veca=zeros(1,dimension);
            vecb=zeros(1,dimension);
            veca(j)=1;
            vecb(j+k)=1;
            f=(j-1)*dimension-(j*(j-1))/2+k;
            jacobian(1,i)=jacobian(1,i)+0.5*ccf(f)*((veca*perp')*(vecb*base')-(vecb*perp')*(veca*base'));
        end
    end
end

a=jacobian;

end

function [A,Aeq]=nonlcon(y,s,n,dimension,lb,ub)
%%%%%%%%% 
%
%This function imposes the nonlinear constraint that all the points forming the gait stay within bounds
%
%Inputs:
%
%y: Fourier series coefficients that describe the gait
%s: System file which contains the connection vector field, CCF's and
%   metric data
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
%dimension: Indicates the number of shape variables of the system
%lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
%ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% 
%%%%%%%%%

% % The first step is to obtain a direct transciption parametrization of the gait from the 
% % fourier series parametrization
y1 = path_from_fourier(y,n,dimension);
y2=y1(:);

%b=length(y2);

% A1 and A2 together impose the constraint that all the points forming the gait stay in bounds
A1=y2+lb;
A2=-y2-ub;

A = [A1;A2];


% Make sure the frequency doesn't get changed from 2*pi
Aeq = y(end,:) - 2*pi;

end


function y = path_from_fourier(f,n,dimension)
% Returns the shape space parametrization of the gait at n points when provided with
% the fourier coefficients f. The gait is returned as a self-closed gait
% (i.e. the first and last rows of the output y are the same point).
% Inputs:
%   f: Fourier coefficients that parametrize the gait.
%   n: Number of points that should compose the gait, less the final
%       self-closed point
%   dimension: Number of shape variables for system

    y = zeros(n+1,dimension);
    % Determine time period based on value of fourier frequency
    w = f(end,1);
    T = 2*pi/w;
    % Create time vector at which to evaluate points of gait
    t = linspace(0,T,n+1);
    % Evaluate the shape-space parametrization of the gait at every time
    % value in t
    for j=1:dimension
        for i=1:1:n+1
            y(i,j)=f(1,j)+f(2,j)*cos(w*t(i))+f(3,j)*sin(w*t(i))+f(4,j)*cos(2*w*t(i))+...
                +f(5,j)*sin(2*w*t(i))+f(6,j)*cos(3*w*t(i))+f(7,j)*sin(3*w*t(i))+...
                +f(8,j)*cos(4*w*t(i))+f(9,j)*sin(4*w*t(i));
        end
    end
end

function stop=outfun(y,optimValues,state,stretch,s,handles)
%%%%%%%%% 
%
%This function plots the current state of the gait on the sysplotter GUI
%after every iteration of the optimizer
%
%%%%%%%%% 

n=100;
dimension=length(y(1,:));

% % The if else statement below deletes gaits 2 iterations after they have been plotted
% if optimValues.iteration>2
%     children=get(gca,'children');
%     delete(children(6:10));
% else
% end

for thisAxes = [1:numel(handles.plot_thumbnails.Children)]
    
    axes(handles.plot_thumbnails.Children(thisAxes));

    % The if else statement below fades the gait plotted during the previous iteration
    if optimValues.iteration>1
        children=get(gca,'children');
        for idx = 1:numel(children)

            if iscell(children(idx).UserData) && strcmp(children(idx).UserData{1},'OptimizeTracer')
                children(idx).UserData = {'OptimizeTracer', children(idx).UserData{2}-1};

                if children(idx).UserData{2} == 0

                    delete(children(idx));

                else

                    children(idx).Color=[0.5 0.5 0.5];
                    children(idx).LineWidth=4;
                end
            end

        end
    %     children(1).Color=[0.5 0.5 0.5];
    %     children(2).Color=[0.5 0.5 0.5];
    %     children(3).Color=[0.5 0.5 0.5];
    %     children(4).Color=[0.5 0.5 0.5];
    %     children(5).Color=[0.5 0.5 0.5];
    % 
    %     children(1).LineWidth=4;
    %     children(2).LineWidth=4;
    %     children(3).LineWidth=4;
    %     children(4).LineWidth=4;
    %     children(5).LineWidth=4;
    else
    end

    % The if else statement below plots the gait after every iteration
    if optimValues.iteration>0
        y1 = path_from_fourier(y,n,dimension);
        hold on
        if stretch
            stretchnames = {'stretch','surface'};
            stretchname = stretchnames{stretch};

            [x_temp,y_temp,z_temp] = s.convert.(stretchname).old_to_new_points(y1(:,1),y1(:,2));
        else
            x_temp = y1(:,1);
            y_temp = y1(:,2);
            z_temp = zeros(size(y1(:,1)));
        end
        handle1=line('XData',x_temp,'YData',y_temp,'ZData',z_temp,'color','k','linewidth',3,'UserData',{'OptimizeTracer',2}); %#ok<NASGU>
        %plot_dir_arrows(y1(:,1),y1(:,2),2,'Color',[0 0 0],'LineWidth',3);
    else
    end

    % % Use this version if you're writing a video and have the system plotting
    % % more frequently in the difffmincon call
    % if optimValues.iteration>0
    %     y1 = path_from_fourier(y,n,dimension);
    %     figure(5);
    %     subplot(1,2,1)
    %     delete(findobj(gca,'Type','Line'));
    %     handle1=plot(y1(:,1),y1(:,2),'k','linewidth',3);
    %     plot_dir_arrows(y1(:,1),y1(:,2),2,'Color',[0 0 0],'LineWidth',3);
    %     xlabel('\alpha_1')
    %     ylabel('\alpha_2')
    %     
    %     subplot(1,2,2)
    %     if optimValues.iteration > 1
    %         fig = gcf;
    %         axObjs = fig.Children;
    %         dataObjs = axObjs(1).Children;
    %         iterations = [dataObjs(1).XData, optimValues.iteration];
    %         fvals = [dataObjs(1).YData, optimValues.fval];
    %     else
    %         iterations = optimValues.iteration;
    %         fvals = optimValues.fval;
    %     end
    %     
    %     plot(iterations,fvals,'bo-')
    %     xlabel('Optimizer Iteration')
    %     ylabel('Efficiency')
    %     title('Efficiency per iteration')
    % end
end

pause(0.05)
stop=false;
end

function [net_disp_orig, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,tspan,ConnectionEval,IntegrationMethod,resolution)
% Evaluate the displacement and cost for the gait specified in the
% structure GAIT when executed by the system described in the structure
% S.
%
% S should be a sysplotter 's' structure loaded from a file
% sysf_FILENAME_calc.mat (note the _calc suffix)
%
% P should be a structure with fields "phi_def" and "dphi_def", returning a
% vector of shapes and shape velocities respectively. If it is not
% convenient to analytically describe the shape velocity function,
% gait.dphi should be defined as 
%
% p.dphi =  @(t) jacobianest(@(T) p.phi (T),t)
%
% as is done automatically by sysplotter, but note that this will be slower
% than specifying a function directly
%
% ConnectionEval can specify whether the local connection should be generated from
% its original function definiton, or by interpolation into the evaluated
% matrix, but is optional. Valid options are 'functional' or
% 'interpolated'. Defaults to "interpolated", which significantly faster
% when calculating the local connection or metric from scratch takes
% appreciable computational time
%
% IntegrationMethod can specify whether ODE45 or a fixed point
% (euler-exponential) integration method should be employed. Defaults to
% ODE, fixed point code is experimental.
%
% RESOLUTION specifies the number of points for fixed-point resolution
% evaluation. A future option may support autoconvergence, but ODE
% performance with interpolated evaluation appears to be fast enough that
% refinement of fixed-point performance is on hold.
	

	% if no ConnectionEval method is specified, default to interpolated
	if ~exist('ConnectionEval','var')
		ConnectionEval = 'interpolated';
	end
    
    % if no IntegrationMethod is specified, default to ODE
	if ~exist('IntegrationMethod','var')
		IntegrationMethod = 'ODE';
	end

    % if no resolution is specified, default to 100 (this only affects
    % fixed_step integration)
	if ~exist('resolution','var')
		resolution = 100;
	end

    
    
	switch IntegrationMethod
		
		case 'fixed_step'
			
			[net_disp_orig, cost] = fixed_step_integrator(s,p,tspan,ConnectionEval,resolution);
        
        case 'ODE'

            % Calculate the system motion over the gait
            sol = ode45(@(t,y) helper_function(t,y,s,p,ConnectionEval),tspan,[0 0 0 0]');

            % Extract the final motion
            disp_and_cost = deval(sol,tspan(end));

            net_disp_orig = disp_and_cost(1:3);
            cost = disp_and_cost(end);
            
        otherwise
			error('Unknown method for integrating motion');
	end

	
	% Convert the final motion into its representation in optimal
	% coordinates
	startshape = p.phi_def(0);
	startshapelist = num2cell(startshape);
	beta_theta = interpn(s.grid.eval{:},s.B_optimized.eval.Beta{3},startshapelist{:},'spline');
	net_disp_opt = [cos(beta_theta) sin(beta_theta) 0;...
		-sin(beta_theta) cos(beta_theta) 0;...
		0 0 1]*net_disp_orig;

	
end

% Evaluate the body velocity and cost velocity (according to system metric)
% at a given time
function [xi, dcost] = get_velocities(t,s,gait,ConnectionEval)

	% Get the shape and shape derivative at the current time
	shape = gait.phi_def(t);
	shapelist = num2cell(shape);
	dshape = gait.dphi_def(t);
    ddshape = gait.ddphi_def(t);
	
	% Get the local connection and metric at the current time, in the new coordinates
	switch ConnectionEval
		case 'functional'
			
			A = s.A_num(shapelist{:})./s.A_den(shapelist{:});
			
            switch s.system_type
                case 'drag'
                    M = s.metric(shapelist{:});
                case 'inertia'
                    error('Functional ConnectionEval method not supported for inertia systems!')
            end

		case 'interpolated'
			
			A = -cellfun(@(C) interpn(s.grid.eval{:},C,shapelist{:},'spline'),...
                s.vecfield.eval.content.Avec);
			
            metric =  cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                shapelist{:},'spline'),s.metricfield.metric_eval.content.metric);
            switch s.costfunction
                %If our cost is from the coordinate space
                case {'pathlength coord','acceleration coord'}
                    %The reimannian metric is the identity
                    metric = eye(size(metric));
                %If our cost is inertial
                case {'torque','covariant acceleration','power quality'}
                    %Calculate mass matrix
                    M_a = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
                        shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
                    %And mass matrix derivative
                    dM_alphadalpha = cell(size(shapelist));
                    for i = 1:length(shapelist)
                        dM_alphadalpha{i} = cellfun(@(C) interpn(s.grid.coriolis_eval{:},C,...
                            shapelist{:},'spline'),s.coriolisfield.coriolis_eval.content.dM_alphadalpha{i});
                    end
            end
			
		otherwise
			error('Unknown method for evaluating local connection');
	end
	
	% Get the body velocity at the current time
	%t;
    xi = - A * dshape(:);

    switch s.costfunction
        case {'pathlength metric','pathlength coord'}
            dcost = sqrt(dshape(:)'*metric*dshape(:));
        case 'pathlength metric2'
            dcost = sqrt(dshape(:)'*metric*metric*dshape(:));
        case 'torque'
            dcost = torque_cost(M_a,dM_alphadalpha,shape,dshape,ddshape,metric);
        case 'covariant acceleration'
            dcost = acceleration_cost(M_a,dM_alphadalpha,shape,dshape,ddshape,metric);
        case 'acceleration coord'
            dcost = ddshape(:)'*metric*ddshape(:);
        case 'power quality'
            dcost = power_quality_cost(M_a,dM_alphadalpha,shape,dshape,ddshape);
    end
	
end

function dcost = torque_cost(M,dM_alphadalpha,shape,dshape,ddshape,metric)
% Calculates the incremental cost for an inertial system where cost is torque squared.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system

    % Start by calculating the coriolis matrix
    C = calc_coriolis_matrix(dM_alphadalpha,shape,dshape);
    % Calculate the torque for this instant of time and return the inner
    % product of the torque with itself
    dtau = M*ddshape(:) + C;
    dcost = dtau'*metric*dtau;
end

function dcost = acceleration_cost(M,dM_alphadalpha,shape,dshape,ddshape,metric)
% Calculates the incremental cost for an inertial system where cost is covariant acceleration.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system
    C = calc_coriolis_matrix(dM_alphadalpha,shape,dshape);
    cov_acc = ddshape(:) + inv(M)*C;
    dcost = cov_acc'*metric*cov_acc;

end

function dcost = power_quality_cost(M,dM_alphadalpha,shape,dshape,ddshape)
% Calculates the incremental cost for an inertial system where cost is power quality.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system

    % Start by calculating the coriolis matrix
    C = calc_coriolis_matrix(dM_alphadalpha,shape,dshape);
    % Calculate the torque for this instant of time 
    dtau = M*ddshape(:) + C;
    % Calculate power quality
    dcost = (dshape(:)'*dtau)^2 - ((dshape(:)').^2*dtau.^2);
    dcost = dcost + 100;
end

% Function to integrate up system velocities using a fixed-step method
function [net_disp_orig, cost] = fixed_step_integrator(s,gait,tspan,ConnectionEval,resolution)

	% Duplicate 'resolution' to 'res' if it is a number, or place res at a
	% starting resolution if an automatic convergence method is selected
	% (automatic convergence not yet enabled)
	default_res = 100;
	if isnumeric(resolution)
		res = resolution;
	elseif ischar(resolution) && strcmp(resolution,'autoconverge')
		res = default_res;
	else
		error('Unexpected value for resolution');
	end
	
	% Generate the fixed points from the time span and resolution
	tpoints = linspace(tspan(1),tspan(2),res);
	tsteps = gradient(tpoints);

	% Evaluate the velocity function at each time
	[xi, dcost] = arrayfun(@(t) get_velocities(t,s,gait,ConnectionEval),tpoints,'UniformOutput',false);
	
	
	%%%%%%%
	% Integrate cost and displacement into final values
	
	%%
	% Exponential integration for body velocity
	
	% Exponentiate each velocity over the corresponding time step
	expXi = cellfun(@(xi,timestep) se2exp(xi*timestep),xi,num2cell(tsteps),'UniformOutput',false);
	
	% Start off with zero position and displacement
	net_disp_matrix = eye(size(expXi{1}));
	
	% Loop over all the time steps from 1 to n-1, multiplying the
	% transformation into the current displacement
	for i = 1:(length(expXi)-1)
		
		net_disp_matrix = net_disp_matrix * expXi{i};
		
	end
	
	% De-matrixafy the result
	g_theta = atan2(net_disp_matrix(2,1),net_disp_matrix(1,1));
	g_xy = net_disp_matrix(1:2,3);
	
	net_disp_orig = [g_xy;g_theta];
	
	%%
	% Trapezoidal integration for cost
	dcost = [dcost{:}];
	cost = trapz(tpoints,dcost);

end


% Function to evaluate velocity and differential cost at each time for ODE
% solver
function dX = helper_function(t,X,s,gait,ConnectionEval)

	% X is the accrued displacement and cost

	[xi, dcost] = get_velocities(t,s,gait,ConnectionEval);
		
	% Rotate body velocity into world frame
	theta = X(3);
	v = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*xi;
		
	% Combine the output
	dX = [v;dcost];
	

end

function expXi = se2exp(xi)

	% Make sure xi is a column
	xi = xi(:);

	% Special case non-rotating motion
	if xi(3) == 0
		
		expXi = [eye(2) xi(1:2); 0 0 1];
		
	else
		
		z_theta = xi(3);
		
		z_xy = 1/z_theta * [sin(z_theta), 1-cos(z_theta); cos(z_theta)-1, sin(z_theta)] * xi(1:2);
		
		expXi = [ [cos(z_theta), -sin(z_theta); sin(z_theta), cos(z_theta)], z_xy;
			0 0 1];
		
	end


end

% function [g_end_orig,g_end_opt, cost_end] = extract_displacement_and_cost(datafile)
% % Extract the displacement and cost data from a sysf_...shchf_....mat file
% 
% % Load the target file
% load(datafile,'p')
% 
% % Prime arrays to hold the net displacement (in original and optimal
% % coordinates) and cost from each shape change in the file. p.G_locus_full is
% % single-level cell array of structures, each of which holds the
% % information for one gait (with all segments concatenated)
% g_end_orig = zeros(numel(p.G_locus_full),3);
% g_end_opt = g_end_orig;
% cost_end = zeros(numel(p.G_locus_full,1)); % If distance metric was not specified, euclidean metric in the parameters was assumed
% 
% % Loop over each shape change
% for i = 1:numel(p.G_locus_full)
% 	
% 	% Extract the final values for the relevant parameters
% 	g_end_orig(i,:) = p.G_locus_full{i}.G(end,:); 
% 	g_end_opt(i,:) = p.G_locus_full{i}.G_opt(end,:); 
% 	cost_end(i) = p.G_locus_full{i}.S(end);
% end
% end

function [grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g)
% Calculates the gradient of the shape position, velocity, and acceleration
% with respect to the fourier coefficients.
% Inputs:
%   n: Number of points composing the gait when using shape-space
%       parametrization
%   y: Set of fourier coefficients parametrizing the gait
%   g: Time period of the gait

% Get the fourier frequency, number of shape variables, and number of
% fourier coefficientsinertia_cost_gradient
w = y(end,:);
dim = size(y,2);
num_coeffs = size(y,1);
% Initialize cell array of function handles to hold the partials of the
% shape variables with respect to the fourier coefficients
grad_alpha = cell(num_coeffs,dim);
grad_alphadot = cell(num_coeffs,dim);
grad_alphaddot = cell(num_coeffs,dim);

for i = 1:num_coeffs 
    for j = 1:dim
        % Coefficient a_0 is a lone scalar, so partials with respect to it
        % are zero for the dotted terms and 1 for alpha
        if i == 1
            grad_alpha{i,j} = @(t) [0*(1:j-1), 1, 0*(j+1:dim)];
            grad_alphadot{i,j} = @(t) zeros(1,dim);
            grad_alphaddot{i,j} = @(t) zeros(1,dim);
            continue
        elseif i == num_coeffs % Partial w.r.t. frequency
            grad_alpha{i,j} = @(t) [0*(1:j-1), t*(-y(2,j)*sin(w(j)*t) + y(3,j)*cos(w(j)*t) - ...
                                               2*y(4,j)*sin(2*w(j)*t) + 2*y(5,j)*cos(2*w(j)*t) - ...
                                               3*y(6,j)*sin(3*w(j)*t) + 3*y(7,j)*cos(3*w(j)*t) - ...
                                               4*y(8,j)*sin(4*w(j)*t) + 4*y(9,j)*cos(4*w(j)*t)), ...
                                    0*(j+1:dim)];
            
            grad_alphadot{i,j} = @(t) [0*(1:j-1), -y(2,j)*(w(j)*t*cos(w(j)*t)+sin(w(j)*t)) + y(3,j)*(-w(j)*t*sin(w(j)*t)+cos(w(j)*t)) + ...
                                                  -2*y(4,j)*(2*w(j)*t*cos(2*w(j)*t)+sin(2*w(j)*t)) + 2*y(5,j)*(-2*w(j)*t*sin(2*w(j)*t)+cos(2*w(j)*t)) + ...
                                                  -3*y(6,j)*(3*w(j)*t*cos(3*w(j)*t)+sin(3*w(j)*t)) + 3*y(7,j)*(-3*w(j)*t*sin(3*w(j)*t)+cos(3*w(j)*t)) + ...
                                                  -4*y(8,j)*(4*w(j)*t*cos(4*w(j)*t)+sin(4*w(j)*t)) + 4*y(9,j)*(-4*w(j)*t*sin(4*w(j)*t)+cos(4*w(j)*t)), ...
                                      0*(j+1:dim)];
            
            grad_alphaddot{i,j} = @(t) [0*(1:j-1), -y(2,j)*w(j)*(-t*w(j)*sin(w(j)*t)+2*cos(w(j)*t)) - y(3,j)*w(j)*(t*w(j)*cos(w(j)*t)+2*sin(w(j)*t)) + ...
                                                   -4*y(4,j)*w(j)*(-2*t*w(j)*sin(2*w(j)*t)+2*cos(2*w(j)*t)) - 4*y(5,j)*w(j)*(2*t*w(j)*cos(2*w(j)*t)+2*sin(2*w(j)*t)) + ...
                                                   -9*y(6,j)*w(j)*(-3*t*w(j)*sin(3*w(j)*t)+2*cos(3*w(j)*t)) - 9*y(7,j)*w(j)*(3*t*w(j)*cos(3*w(j)*t)+2*sin(3*w(j)*t)) + ...
                                                   -16*y(8,j)*w(j)*(-4*t*w(j)*sin(4*w(j)*t)+2*cos(4*w(j)*t)) - 16*y(9,j)*w(j)*(4*t*w(j)*cos(4*w(j)*t)+2*sin(4*w(j)*t)), ...
                                       0*(j+1:dim)];
            continue
        end
        % For partial alpha, a_n is associated with cosine and b_n is
        % associated with sine; a_n terms are every second row entry in y
        % with the b_n terms in between
        if mod(i,2) == 0
            trig = @cos;
        else
            trig = @sin;
        end
        % mult comes from the multiplier of the natural frequency for
        % increasing fourier coefficients
        mult = floor(i/2);
        
        grad_alpha{i,j} = @(t) [0*(1:j-1), trig(mult*w(j)*t), 0*(j+1:dim)];
        % For partial alphadot, a_n is associated with sine and b_n is
        % associated with cosine; the a_n terms are every second row entry 
        % in y with the b_n terms in between
        if mod(i,2) == 0
            trig = @sin;
        else
            trig = @cos;
        end
        
        grad_alphadot{i,j} = @(t) [0*(1:j-1), (-1)^(i-1)*mult*w(j)*trig(mult*w(j)*t), 0*(j+1:dim)];
        
        % For partial alphaddot, a_n is associated with cosine and b_n is
        % associated with sine
        if mod(i,2) == 0
            trig = @cos;
        else
            trig = @sin;
        end
        
        grad_alphaddot{i,j} = @(t) [0*(1:j-1), -mult^2*w(j)^2*trig(mult*w(j)*t), 0*(j+1:dim)];
    end
end
end

function cost_grad = inertia_cost_gradient(s,n,y,g,p,EvaluationMethod)
% Calculates the gradient of cost for inertial systems.
% Inputs:
%   s: System struct used by sysplotter.
%   n: Number of points at which gait should be evaluated in the shape
%       space.
%   y: Fourier coefficients that parametrize the gait.
%   g: Time period over which gait is executed.
%   p: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   EvaluationMethod: String representing how the gradient of cost should
%   be calculated; provide 'discrete' to evaluate at 100 discrete values or
%   'ode45' if you would like the gradient of cost to be integrated using
%   ode45. Other values will result in error.

    % Contribution to gradient from the movement of each point due to
    % change in fourier coefficients
    [grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g);

    cost_grad = zeros(size(grad_alpha));
    tspan = [0 g];
    if strcmpi(EvaluationMethod,'discrete')
        num_pts = 100;
        t_pts = linspace(0,g,num_pts);
        del_t = t_pts(2) - t_pts(1);
        switch s.costfunction
            case 'torque'
                for k = 1:length(t_pts)
                    del_cost = inertia_gradient_helper(t_pts(k),[],s,p,grad_alpha,grad_alphadot,grad_alphaddot);
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'covariant acceleration'
                for k = 1:length(t_pts)
                    del_cost = acceleration_gradient_helper(t_pts(k),[],s,p,grad_alpha,grad_alphadot,grad_alphaddot);
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'acceleration coord'
                for k = 1:length(t_pts)
                    del_cost = accelerationcoord_gradient_helper(t_pts(k),[],s,p,grad_alphaddot);
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'power quality'
                for k = 1:length(t_pts)
                    del_cost = powerquality_gradient_helper(t_pts(k),[],s,p,grad_alpha,grad_alphadot,grad_alphaddot);
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
        end
        % Reset gradient of fourier frequency to be zero to prevent changes
        % to it
        cost_grad(end,:) = 0;
    elseif strcmpi(EvaluationMethod,'ode45')
        switch s.costfunction
            case 'torque'
                sol = ode45(@(t,y) inertia_gradient_helper(t,y,s,p,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
            case 'covariant acceleration'
                sol = ode45(@(t,y) acceleration_gradient_helper(t,y,s,p,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
            case 'acceleration coord'
                sol = ode45(@(t,y) accelerationcoord_gradient_helper(t,y,s.p,grad_alphaddot),tspan,cost_grad);
            case 'power quality'
                sol = ode45(@(t,y) powerquality_gradient_helper(t,y,s,p,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
        end
        % Extract the final motion
        cost_grad = reshape(deval(sol,tspan(end)),size(cost_grad));
        % Reset gradient of fourier frequency to be zero to prevent changes
        % to it
        cost_grad(end,:) = 0;
    else
        error('Untenable option provided for EvaluationMethod!')
    end
end

function validate_shape_gradient(n,y,g,grad_alphaddot,grad_alphadot,grad_alpha) %#ok<DEFNU>
% Function that helps validate that the gradient of shape position,
% velocity, and acceleration are correctly calculated. Difference between
% the input gradients and calculation-verified gradients are printed to the
% terminal. Should be very close to zero.
% Inputs:
%   n: Number of points at which gait should be evaluated in shape space
%   y: Fourier coefficients that parametrize the gait
%   g: Time period over which gait is executed
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

    % Select 10 random times at which to evaluate the gradient of shape
    t = g*rand(1,10);
    w1 = y(end,1); % Frequency of Fourier transform
    w2 = y(end,2);
    
    % Evaluate the gradient at each random time against a hard-coded
    % gradient calculation
    for i = 1:length(t)
        grad_alpha_eval = cellfun(@(C) C(t(i)), grad_alpha, 'UniformOutput', false);
        grad_alpha1_eval = cell2mat(grad_alpha_eval(:,1));
        grad_alpha2_eval = cell2mat(grad_alpha_eval(:,2));
        grad_alpha1_calc = [1,0;cos(w1*t(i)),0;sin(w1*t(i)),0;cos(2*w1*t(i)),0;sin(2*w1*t(i)),0;cos(3*w1*t(i)),0;sin(3*w1*t(i)),0;cos(4*w1*t(i)),0;sin(4*w1*t(i)),0;0,0];
        grad_alpha2_calc = [0,1;0,cos(w2*t(i));0,sin(w2*t(i));0,cos(2*w2*t(i));0,sin(2*w2*t(i));0,cos(3*w2*t(i));0,sin(3*w2*t(i));0,cos(4*w2*t(i));0,sin(4*w2*t(i));0,0];
        grad_alpha1_err = grad_alpha1_eval - grad_alpha1_calc
        grad_alpha2_err = grad_alpha2_eval - grad_alpha2_calc
        
        grad_alphadot_eval = cellfun(@(C) C(t(i)), grad_alphadot, 'UniformOutput', false);
        grad_alphadot1_eval = cell2mat(grad_alphadot_eval(:,1));
        grad_alphadot2_eval = cell2mat(grad_alphadot_eval(:,2));
        grad_alphadot1_calc = [0,0;-w1*sin(w1*t(i)),0;w1*cos(w1*t(i)),0;-2*w1*sin(2*w1*t(i)),0;2*w1*cos(2*w1*t(i)),0;-3*w1*sin(3*w1*t(i)),0;3*w1*cos(3*w1*t(i)),0;-4*w1*sin(4*w1*t(i)),0;4*w1*cos(4*w1*t(i)),0;0,0];
        grad_alphadot2_calc = [0,0;0,-w2*sin(w2*t(i));0,w2*cos(w2*t(i));0,-2*w2*sin(2*w2*t(i));0,2*w2*cos(2*w2*t(i));0,-3*w2*sin(3*w2*t(i));0,3*w2*cos(3*w2*t(i));0,-4*w2*sin(4*w2*t(i));0,4*w2*cos(4*w2*t(i));0,0];
        grad_alphadot1_err = grad_alphadot1_eval - grad_alphadot1_calc;
        grad_alphadot2_err = grad_alphadot2_eval - grad_alphadot2_calc;
        
        grad_alphaddot_eval = cellfun(@(C) C(t(i)), grad_alphaddot, 'UniformOutput', false);
        grad_alphaddot1_eval = cell2mat(grad_alphaddot_eval(:,1));
        grad_alphaddot2_eval = cell2mat(grad_alphaddot_eval(:,2));
        grad_alphaddot1_calc = [0,0;-w1^2*cos(w1*t(i)),0;-w1^2*sin(w1*t(i)),0;-4*w1^2*cos(2*w1*t(i)),0;-4*w1^2*sin(2*w1*t(i)),0;-9*w1^2*cos(3*w1*t(i)),0;-9*w1^2*sin(3*w1*t(i)),0;-16*w1^2*cos(4*w1*t(i)),0;-16*w1^2*sin(4*w1*t(i)),0;0,0];
        grad_alphaddot2_calc = [0,0;0,-w2^2*cos(w2*t(i));0,-w2^2*sin(w2*t(i));0,-4*w2^2*cos(2*w2*t(i));0,-4*w2^2*sin(2*w2*t(i));0,-9*w2^2*cos(3*w2*t(i));0,-9*w2^2*sin(3*w2*t(i));0,-16*w2^2*cos(4*w2*t(i));0,-16*w2^2*sin(4*w2*t(i));0,0];
        grad_alphaddot1_err = grad_alphaddot1_eval - grad_alphaddot1_calc;
        grad_alphaddot2_err = grad_alphaddot2_eval - grad_alphaddot2_calc;
    end
end

function validate_cost_gradient(s,n,y,g,p)
% Function to help verify that the cost gradient calculated by
% inertia_gradient_helper is valid. Values returned by function and
% individually calculated costs are printed to the terminal along with the
% difference between the two methods. Should be relatively small for full
% range of fourier coefficients.
% Inputs:
%   s: System struct used by sysplotter.
%   n: Number of points at which gait should be evaluated in the shape
%       space.
%   y: Fourier coefficients that parametrize the gait.
%   g: Time period over which gait is executed.
%   p: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t

% Set the delta in the fourier coefficients between individual cost
% evaluations
fourier_delta = 0.001;
[grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g);
% Go through each of the fourier coefficients and test how applying
% fourier_delta to each of them matches with the cost gradient calculated
% by inertia_gradient_helper
for fourier_test = 1:numel(y)
    % Perturb the fourier coefficient to be tested by fourier_delta and
    % create a parametrization for calculating the point in the gait at a
    % time t
    y2 = y;
    y2(fourier_test) = y2(fourier_test) + fourier_delta;
    w1 = y2(end,1); % Frequency of Fourier transform
    w2 = y2(end,2);
    p2.phi_def = @(t) [y2(1,1)+y2(2,1)*cos(w1*t)+y2(3,1)*sin(w1*t)+y2(4,1)*cos(2*w1*t)+...
                +y2(5,1)*sin(2*w1*t)+y2(6,1)*cos(3*w1*t)+y2(7,1)*sin(3*w1*t)+...
                +y2(8,1)*cos(4*w1*t)+y2(9,1)*sin(4*w1*t),y2(1,2)+y2(2,2)*cos(w2*t)+y2(3,2)*sin(w2*t)+y2(4,2)*cos(2*w2*t)+...
                +y2(5,2)*sin(2*w2*t)+y2(6,2)*cos(3*w2*t)+y2(7,2)*sin(3*w2*t)+...
                +y2(8,2)*cos(4*w2*t)+y2(9,2)*sin(4*w2*t)]; % function parametrizing the gait as a function of time
    p2.dphi_def = @(t) [-(w1)*y2(2,1)*sin(w1*t)+(w1)*y2(3,1)*cos(w1*t)-2*(w1)*y2(4,1)*sin(2*w1*t)+...
                +2*(w1)*y2(5,1)*cos(2*w1*t)-3*(w1)*y2(6,1)*sin(3*w1*t)+3*(w1)*y2(7,1)*cos(3*w1*t)+...
                -4*(w1)*y2(8,1)*sin(4*w1*t)+4*(w1)*y2(9,1)*cos(4*w1*t),-(w2)*y2(2,2)*sin(w2*t)+...
                (w2)*y2(3,2)*cos(w2*t)-2*(w2)*y2(4,2)*sin(2*w2*t)+...
                +2*(w2)*y2(5,2)*cos(2*w2*t)-3*(w2)*y2(6,2)*sin(3*w2*t)+3*(w2)*y2(7,2)*cos(3*w2*t)+...
                -4*(w2)*y2(8,2)*sin(4*w2*t)+4*(w2)*y2(9,2)*cos(4*w2*t)]; % Shape space velocity as a function of time
    p2.ddphi_def = @(t) [-(w1)^2*y2(2,1)*cos(w1*t)-(w1)^2*y2(3,1)*sin(w1*t)-4*(w1)^2*y2(4,1)*cos(2*w1*t)+...
                -4*(w1)^2*y2(5,1)*sin(2*w1*t)-9*(w1)^2*y2(6,1)*cos(3*w1*t)-9*(w1)^2*y2(7,1)*sin(3*w1*t)+...
                -16*(w1)^2*y2(8,1)*cos(4*w1*t)-16*(w1)^2*y2(9,1)*sin(4*w1*t),-(w2)^2*y2(2,2)*cos(w2*t)+...
                -(w2)^2*y2(3,2)*sin(w2*t)-4*(w2)^2*y2(4,2)*cos(2*w2*t)+...
                -4*(w2)^2*y2(5,2)*sin(2*w2*t)-9*(w2)^2*y2(6,2)*cos(3*w2*t)-9*(w2)^2*y2(7,2)*sin(3*w2*t)+...
                -16*(w2)^2*y2(8,2)*cos(4*w2*t)-16*(w2)^2*y2(9,2)*sin(4*w2*t)]; % Shape space accel. as a function of time
            
    % Perturb again in the opposite direction
    y2 = y;
    y2(fourier_test) = y2(fourier_test) - fourier_delta;
    w1 = y2(end,1); % Frequency of Fourier transform
    w2 = y2(end,2);
    p3.phi_def = @(t) [y2(1,1)+y2(2,1)*cos(w1*t)+y2(3,1)*sin(w1*t)+y2(4,1)*cos(2*w1*t)+...
                +y2(5,1)*sin(2*w1*t)+y2(6,1)*cos(3*w1*t)+y2(7,1)*sin(3*w1*t)+...
                +y2(8,1)*cos(4*w1*t)+y2(9,1)*sin(4*w1*t),y2(1,2)+y2(2,2)*cos(w2*t)+y2(3,2)*sin(w2*t)+y2(4,2)*cos(2*w2*t)+...
                +y2(5,2)*sin(2*w2*t)+y2(6,2)*cos(3*w2*t)+y2(7,2)*sin(3*w2*t)+...
                +y2(8,2)*cos(4*w2*t)+y2(9,2)*sin(4*w2*t)]; % function parametrizing the gait as a function of time
    p3.dphi_def = @(t) [-(w1)*y2(2,1)*sin(w1*t)+(w1)*y2(3,1)*cos(w1*t)-2*(w1)*y2(4,1)*sin(2*w1*t)+...
                +2*(w1)*y2(5,1)*cos(2*w1*t)-3*(w1)*y2(6,1)*sin(3*w1*t)+3*(w1)*y2(7,1)*cos(3*w1*t)+...
                -4*(w1)*y2(8,1)*sin(4*w1*t)+4*(w1)*y2(9,1)*cos(4*w1*t),-(w2)*y2(2,2)*sin(w2*t)+...
                (w2)*y2(3,2)*cos(w2*t)-2*(w2)*y2(4,2)*sin(2*w2*t)+...
                +2*(w2)*y2(5,2)*cos(2*w2*t)-3*(w2)*y2(6,2)*sin(3*w2*t)+3*(w2)*y2(7,2)*cos(3*w2*t)+...
                -4*(w2)*y2(8,2)*sin(4*w2*t)+4*(w2)*y2(9,2)*cos(4*w2*t)]; % Shape space velocity as a function of time
    p3.ddphi_def = @(t) [-(w1)^2*y2(2,1)*cos(w1*t)-(w1)^2*y2(3,1)*sin(w1*t)-4*(w1)^2*y2(4,1)*cos(2*w1*t)+...
                -4*(w1)^2*y2(5,1)*sin(2*w1*t)-9*(w1)^2*y2(6,1)*cos(3*w1*t)-9*(w1)^2*y2(7,1)*sin(3*w1*t)+...
                -16*(w1)^2*y2(8,1)*cos(4*w1*t)-16*(w1)^2*y2(9,1)*sin(4*w1*t),-(w2)^2*y2(2,2)*cos(w2*t)+...
                -(w2)^2*y2(3,2)*sin(w2*t)-4*(w2)^2*y2(4,2)*cos(2*w2*t)+...
                -4*(w2)^2*y2(5,2)*sin(2*w2*t)-9*(w2)^2*y2(6,2)*cos(3*w2*t)-9*(w2)^2*y2(7,2)*sin(3*w2*t)+...
                -16*(w2)^2*y2(8,2)*cos(4*w2*t)-16*(w2)^2*y2(9,2)*sin(4*w2*t)]; % Shape space accel. as a function of time
     
    % Get the shape and shape derivative at a random time for each gait
    t = g*rand(1);
	shape = p3.phi_def(t);
	shapelist = num2cell(shape);
	dshape = p3.dphi_def(t);
    ddshape = p3.ddphi_def(t);
    shape_delta = p2.phi_def(t);
	shapelist_delta = num2cell(shape_delta);
	dshape_delta = p2.dphi_def(t);
    ddshape_delta = p2.ddphi_def(t);
    
    % Get mass matrices at both locations
    M = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    M_delta = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist_delta{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    
    % Get partial mass matrices at both locations
    dM_alphadalpha = calc_partial_mass(s,shapelist);
    dM_alphadalpha_delta = calc_partial_mass(s,shapelist_delta);
    
    % Calculate the cost using the two different gaits
    cost = torque_cost(M,dM_alphadalpha,shape,dshape,ddshape);
    cost_delta = torque_cost(M_delta,dM_alphadalpha_delta,shape_delta,dshape_delta,ddshape_delta);
    % Calculate what the gradient of cost is for this particular point
    cost_grad = inertia_gradient_helper(t,[],s,p,grad_alpha,grad_alphadot,grad_alphaddot);
    cost_grad = reshape(cost_grad,size(y));
    cost_grad_rel = cost_grad(fourier_test)
    cost_grad_calc = (cost_delta-cost)/(2*fourier_delta)
    
    % Find what the difference is between cost_grad and the costs evaluated
    % at distance fourier_delta
    err = cost_grad_rel - cost_grad_calc
end

end

% This function could be helpful to verify individual components of the
% torque calculation if things aren't working as expected
% function inertia_gradient_validation(s,n,y,g,p,p_delta,fourier_test,fourier_delta,alpha_test)
%     [grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g);
%     
%     % Get the shape and shape derivative at a random time for each gait
%     t = g*rand(1);
% 	shape = p.phi_def(t);
% 	shapelist = num2cell(shape);
% 	dshape = p.dphi_def(t);
%     ddshape = p.ddphi_def(t);
%     shape_delta = p_delta.phi_def(t);
% 	shapelist_delta = num2cell(shape_delta);
% 	dshape_delta = p_delta.dphi_def(t);
%     ddshape_delta = p_delta.ddphi_def(t);
%     grad_alpha_eval = cellfun(@(C) C(t), grad_alpha, 'UniformOutput', false);
%     grad_alphadot_eval = cellfun(@(C) C(t), grad_alphadot, 'UniformOutput', false);
%     grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
%     
%     %----- dMdf validation -----%
%     % Get mass matrices at both locations
%     M = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
%         shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
%     M_delta = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
%         shapelist_delta{:},'spline'),s.massfield.mass_eval.content.M_alpha);
%     
%     % Get partial mass at both locations
% %     dM_alphadalpha = cell(size(shapelist));
% %     dM_alphadalpha_delta = cell(size(shapelist));
% %     for i = 1:length(shapelist)
% %         dM_alphadalpha{i} = cellfun(@(C) interpn(s.grid.coriolis_eval{:},C,...
% %             shapelist{:},'spline'),s.coriolisfield.coriolis_eval.content.dM_alphadalpha{i});
% %         dM_alphadalpha_delta{i} = cellfun(@(C) interpn(s.grid.coriolis_eval{:},C,...
% %             shapelist_delta{:},'spline'),s.coriolisfield.coriolis_eval.content.dM_alphadalpha{i});
% %     end
%     dM_alphadalpha = calc_partial_mass(s,shapelist);
%     dM_alphadalpha_delta = calc_partial_mass(s,shapelist_delta);
%     
%     % Extract the partial shape derivative wrt fourier coefficient for the
%     % given fourier coefficient index fourier_test corresponding to the
%     % changed shape variable alpha_test
%     grad_alpha_delta = grad_alpha_eval{fourier_test,alpha_test};
%     dMdf = zeros(size(M));
%     for i = 1:length(dM_alphadalpha)
%         dMdf = dMdf + dM_alphadalpha{i}*grad_alpha_delta(i);
%     end
%     
%     dMdf
%     dMdf_calc = (M_delta-M)/fourier_delta
%     % NOTE: This calculation appears to be correct, or at least close
%     % enough. Accurate to two decimals or so of precision.
%     
%     % ----- d(M*alphaddot)/df validation ------ %
%     % Get M*alphaddot at both locations
%     grad_alphaddot_delta = grad_alphaddot_eval{fourier_test,alpha_test};
%     M_alphaddot = M*ddshape(:);
%     M_alphaddot_delta = M_delta*ddshape_delta(:);
%     dMalphaddot_df = dMdf*ddshape(:) + M*grad_alphaddot_delta(:);
%     
%     dMalphaddot_df
%     dMalphaddot_df_calc = (M_alphaddot_delta - M_alphaddot)/fourier_delta
%     % NOTE: This calculation also appears to be close enough. Accurate to
%     % two decimals of precision.
%     
%     % ----- d/df(dM/dq*qdot) validation ----- %
%     grad_alphadot_delta = grad_alphadot_eval{fourier_test,alpha_test};
%     % Get ddM_alphaddalpha at original location
% %     ddM_alphadalpha = cell(length(shapelist));
% %     for i = 1:numel(ddM_alphadalpha)
% %         ddM_alphadalpha{i} = cellfun(@(C) interpn(s.grid.coriolis_eval{:},C,...
% %             shapelist{:},'spline'),s.coriolisfield.coriolis_eval.content.ddM_alphadalpha{i});
% %     end
%     ddM_alphadalpha = calc_second_partial_mass(s,shapelist);
%     
%     % Get dM_alphadalpha*qdot at both locations
%     C_temp = zeros(length(shapelist));
%     C_temp_delta = zeros(length(shapelist));
%     for i = 1:length(dM_alphadalpha)
%         C_temp = C_temp + dM_alphadalpha{i}*dshape(i);
%         C_temp_delta = C_temp_delta + dM_alphadalpha_delta{i}*dshape_delta(i);
%     end
%     
%     % Now see if that matches with what the gradient calculation says
%     C1_partialgrad = zeros(length(shapelist));
%     C1_shapegrad = zeros(length(shapelist));
%     C1_outergrad = zeros(length(shapelist));
%     del_dM_alphadalpha = cell(size(shapelist));
%     for j = 1:length(shapelist)
%         Cj_temp = zeros(length(shapelist));
%         for k = 1:length(shapelist)
%             Cj_temp = Cj_temp + ddM_alphadalpha{j,k}*grad_alpha_delta(k);
%         end
%         % del_dM_alphadalpha is the result of the interior summation
%         del_dM_alphadalpha{j} = Cj_temp;
%         % C1_partialgrad is the exterior summation
%         C1_partialgrad = C1_partialgrad + Cj_temp*dshape(j);
%         C1_shapegrad = C1_shapegrad + dM_alphadalpha{j}*grad_alphadot_delta(j);
%         C1_outergrad = C1_outergrad + dM_alphadalpha{j}*dshape(j);
%     end
%     C1_grad = (C1_partialgrad + C1_shapegrad)
%     C1_grad_calc = (C_temp_delta - C_temp)/fourier_delta
% 
% end

function del_cost = inertia_gradient_helper(t,X,s,gait,grad_alpha,grad_alphadot,grad_alphaddot)
% Helper function to calculate the gradient of cost for inertia systems.
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alpha_eval = cellfun(@(C) C(t), grad_alpha, 'UniformOutput', false);
    grad_alphadot_eval = cellfun(@(C) C(t), grad_alphadot, 'UniformOutput', false);
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alpha_eval));
    % Get the shape and shape derivative at the current time
	shape = gait.phi_def(t);
	shapelist = num2cell(shape);
	dshape = gait.dphi_def(t);
    ddshape = gait.ddphi_def(t);
    
    [metric,metricgrad] = getMetricGrad(s,shape,grad_alpha_eval);
    
    % Get mass and partial mass matrices
    M = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    dM_alphadalpha = calc_partial_mass(s,shapelist);
    ddM_alphadalpha = calc_second_partial_mass(s,shapelist);
    
    % Regular torque calculation
    C = calc_coriolis_matrix(dM_alphadalpha,shape,dshape);
    tau = M*ddshape(:) + C;
    
    for i = 1:numel(grad_alpha_eval)
        % Partial of shape variables with respect to fourier coefficient i
        del_shape = grad_alpha_eval{i};
        del_dshape = grad_alphadot_eval{i};
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of torque calculation
        % Start with effect of gradient on M_alpha*alphaddot
        M_temp = zeros(length(shapelist));
        for j = 1:length(shapelist)
            M_temp = M_temp + dM_alphadalpha{j}*del_shape(j);
        end
        % Catching for debugging
        try
            M_grad = M_temp*ddshape(:) + M*del_ddshape(:);
        catch
            M_temp
        end
        % Effect of gradient on dM_alphadalpha*alphadot*alphadot
        C1_partialgrad = zeros(length(shapelist));
        C1_shapegrad = zeros(length(shapelist));
        C1_outergrad = zeros(length(shapelist));
        del_dM_alphadalpha = cell(size(shapelist));
        for j = 1:length(shapelist)
            Cj_temp = zeros(length(shapelist));
            for k = 1:length(shapelist)
                Cj_temp = Cj_temp + ddM_alphadalpha{j,k}*del_shape(k);
            end
            del_dM_alphadalpha{j} = Cj_temp;
            C1_partialgrad = C1_partialgrad + Cj_temp*dshape(j);
            C1_shapegrad = C1_shapegrad + dM_alphadalpha{j}*del_dshape(j);
            C1_outergrad = C1_outergrad + dM_alphadalpha{j}*dshape(j);
        end
        C1_grad = (C1_partialgrad + C1_shapegrad)*dshape(:) + ...
            C1_outergrad*del_dshape(:);

        % Effect of gradient on -(1/2)*alphadot'*dM_alphadalpha*alphadot
        C2_grad = zeros(size(shapelist(:)));
        for j = 1:length(shapelist)
            C2_grad(j) = del_dshape(:)'*dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*del_dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*dM_alphadalpha{j}*del_dshape(:);
        end
        
        % Gradient of torque
        del_tau = M_grad + C1_grad - (1/2)*C2_grad;
        del_cost(i) = del_tau(:)'*metric*tau(:)...
                    + tau(:)'*metricgrad{i}*tau(:)...
                    + tau(:)'*metric*del_tau(:);
    end
    del_cost = del_cost(:);
end

function del_cost = acceleration_gradient_helper(t,X,s,gait,grad_alpha,grad_alphadot,grad_alphaddot)
% Helper function to calculate the gradient of covariant acceleration cost
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alpha_eval = cellfun(@(C) C(t), grad_alpha, 'UniformOutput', false);
    grad_alphadot_eval = cellfun(@(C) C(t), grad_alphadot, 'UniformOutput', false);
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alpha_eval));
    % Get the shape and shape derivative at the current time
	shape = gait.phi_def(t);
	shapelist = num2cell(shape);
	dshape = gait.dphi_def(t);
    ddshape = gait.ddphi_def(t);
    
    [metric,metricgrad] = getMetricGrad(s,shape,grad_alpha_eval);
    
    % Get mass and partial mass matrices
    M = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    Malpha_inv = inv(M);
    dM_alphadalpha = calc_partial_mass(s,shapelist);
    ddM_alphadalpha = calc_second_partial_mass(s,shapelist);
    
    % Regular torque calculation
    C = calc_coriolis_matrix(dM_alphadalpha,shape,dshape);
    tau = M*ddshape(:) + C;
    
    %Covariant acceleration calculation
    cov_acc = Malpha_inv*tau;
    
    for i = 1:numel(grad_alpha_eval)
        % Partial of shape variables with respect to fourier coefficient i
        del_shape = grad_alpha_eval{i};
        del_dshape = grad_alphadot_eval{i};
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of torque calculation
        % Start with effect of gradient on M_alpha*alphaddot
        M_temp = zeros(length(shapelist));
        for j = 1:length(shapelist)
            M_temp = M_temp + dM_alphadalpha{j}*del_shape(j);
        end
        
        % Catching for debugging
        try
            M_grad = M_temp*ddshape(:) + M*del_ddshape(:);
        catch
            M_temp
        end
        
        %Gradient of inverse of mass for covariant acc calculation
        Minv_grad = zeros(length(shapelist));
        for j = 1:length(shapelist)
            %Formula from matrix cookbook
            %d(M^-1)=-M^-1*dM*M^-1
            dM_alphainv_dalpha = -Malpha_inv*dM_alphadalpha{j}*Malpha_inv;
            Minv_grad = Minv_grad + dM_alphainv_dalpha*del_shape(j);
        end
        
        % Effect of gradient on dM_alphadalpha*alphadot*alphadot
        C1_partialgrad = zeros(length(shapelist));
        C1_shapegrad = zeros(length(shapelist));
        C1_outergrad = zeros(length(shapelist));
        del_dM_alphadalpha = cell(size(shapelist));
        for j = 1:length(shapelist)
            Cj_temp = zeros(length(shapelist));
            for k = 1:length(shapelist)
                Cj_temp = Cj_temp + ddM_alphadalpha{j,k}*del_shape(k);
            end
            del_dM_alphadalpha{j} = Cj_temp;
            C1_partialgrad = C1_partialgrad + Cj_temp*dshape(j);
            C1_shapegrad = C1_shapegrad + dM_alphadalpha{j}*del_dshape(j);
            C1_outergrad = C1_outergrad + dM_alphadalpha{j}*dshape(j);
        end
        C1_grad = (C1_partialgrad + C1_shapegrad)*dshape(:) + ...
            C1_outergrad*del_dshape(:);

        % Effect of gradient on -(1/2)*alphadot'*dM_alphadalpha*alphadot
        C2_grad = zeros(size(shapelist(:)));
        for j = 1:length(shapelist)
            C2_grad(j) = del_dshape(:)'*dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*del_dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*dM_alphadalpha{j}*del_dshape(:);
        end
        % Gradient of torque
        del_tau = M_grad + C1_grad - (1/2)*C2_grad;
        % Gradient of covariant acc
        del_cov_acc = Minv_grad*tau(:)+Malpha_inv*del_tau(:);
        del_cost(i) = del_cov_acc'*metric*cov_acc...
                    + cov_acc'*metricgrad{i}*cov_acc...
                    + cov_acc'*metric*del_cov_acc;
    end
    del_cost = del_cost(:);
end

function del_cost = accelerationcoord_gradient_helper(t,X,s,gait,grad_alphaddot)
% Helper function to calculate the gradient of shape-space acceleration
% cost
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alphaddot_eval));
    % Get the shape derivative at the current time
    ddshape = gait.ddphi_def(t);
   
    % Regular cost calculation
    cost = sqrt(ddshape(:)'*ddshape(:));
    
    for i = 1:numel(grad_alphaddot_eval)
        % Partial of shape acceleration with respect to fourier coefficient i
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of shape-space acceleration cost
        del_cost(i) = del_ddshape(:)'*ddshape(:)+ddshape(:)'*del_ddshape(:);
    end
    del_cost = del_cost(:);
end

function del_cost = powerquality_gradient_helper(t,X,s,gait,grad_alpha,grad_alphadot,grad_alphaddot)
% Helper function to calculate the gradient of cost for inertia systems.
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alpha_eval = cellfun(@(C) C(t), grad_alpha, 'UniformOutput', false);
    grad_alphadot_eval = cellfun(@(C) C(t), grad_alphadot, 'UniformOutput', false);
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alpha_eval));
    % Get the shape and shape derivative at the current time
	shape = gait.phi_def(t);
	shapelist = num2cell(shape);
	dshape = gait.dphi_def(t);
    ddshape = gait.ddphi_def(t);
    
    % Get mass and partial mass matrices
    M = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    dM_alphadalpha = calc_partial_mass(s,shapelist);
    ddM_alphadalpha = calc_second_partial_mass(s,shapelist);
    
    % Regular torque calculation
    C = calc_coriolis_matrix(dM_alphadalpha,shape,dshape);
    tau = M*ddshape(:) + C;
    
    for i = 1:numel(grad_alpha_eval)
        % Partial of shape variables with respect to fourier coefficient i
        del_shape = grad_alpha_eval{i};
        del_dshape = grad_alphadot_eval{i};
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of torque calculation
        % Start with effect of gradient on M_alpha*alphaddot
        M_temp = zeros(length(shapelist));
        for j = 1:length(shapelist)
            M_temp = M_temp + dM_alphadalpha{j}*del_shape(j);
        end
        % Catching for debugging
        try
            M_grad = M_temp*ddshape(:) + M*del_ddshape(:);
        catch
            M_temp
        end
        % Effect of gradient on dM_alphadalpha*alphadot*alphadot
        C1_partialgrad = zeros(length(shapelist));
        C1_shapegrad = zeros(length(shapelist));
        C1_outergrad = zeros(length(shapelist));
        del_dM_alphadalpha = cell(size(shapelist));
        for j = 1:length(shapelist)
            Cj_temp = zeros(length(shapelist));
            for k = 1:length(shapelist)
                Cj_temp = Cj_temp + ddM_alphadalpha{j,k}*del_shape(k);
            end
            del_dM_alphadalpha{j} = Cj_temp;
            C1_partialgrad = C1_partialgrad + Cj_temp*dshape(j);
            C1_shapegrad = C1_shapegrad + dM_alphadalpha{j}*del_dshape(j);
            C1_outergrad = C1_outergrad + dM_alphadalpha{j}*dshape(j);
        end
        C1_grad = (C1_partialgrad + C1_shapegrad)*dshape(:) + ...
            C1_outergrad*del_dshape(:);

        % Effect of gradient on -(1/2)*alphadot'*dM_alphadalpha*alphadot
        C2_grad = zeros(size(shapelist(:)));
        for j = 1:length(shapelist)
            C2_grad(j) = del_dshape(:)'*dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*del_dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*dM_alphadalpha{j}*del_dshape(:);
        end
        
        % Gradient of torque
        del_tau = M_grad + C1_grad - (1/2)*C2_grad;
        % Gradient of (P1+P2)^2, aka gradient of v'*tau*v'*tau where v
        % is the shape velocity
        qual_1 = del_dshape(:)'*tau*dshape(:)'*tau + ...
                 dshape(:)'*del_tau*dshape(:)'*tau + ...
                 dshape(:)'*tau*del_dshape(:)'*tau + ...
                 dshape(:)'*tau*dshape(:)'*del_tau;
        % Gradient of P1^2 + P2^2
        qual_2 = (2*del_dshape(:).*dshape(:))'*(tau.*tau) + ...
                 (dshape(:).*dshape(:))'*(2*del_tau.*tau);
        % Gradient of power quality cost wrt this fourier coeff
        del_cost(i) = qual_1 - qual_2;
    end
    del_cost = del_cost(:);
end

function [metric,metricgrad] = getMetricGrad(s,shape,grad_alpha)
%Returns metric at a shape position, and gradient of metric w.r.t. fourier
%coefficients

%s - structure containing metric function
%shape - array of shape values
%grad_alpha - gradient of shape values w.r.t. fourier coefficients


    %Number of shape variables
    dimension = numel(shape);
    %Step size for small step approximation on metric gradient
    shapestep = 0.0001;
    
    %Declare empty grid for metric interpolation
    interpmetricgrid=cell(1,dimension);
    for j=1:dimension
        interpmetricgrid{j} = s.grid.metric_eval{j,1};
    end
    
    %Calculate metric at shape position
    shapecell = num2cell(shape);
    
    %Initialize metric gradient to zero matrix for all fourier coefficients
    metricgrad = repmat({zeros(dimension)},size(grad_alpha));
    %If cost function is in coordinate space, metric is always identity
    if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'acceleration coord')
        metric = eye(dimension);
        return
    end
    
    metric = zeros(dimension);
    for i = 1:dimension
        for j = 1:dimension
            metric(i,j) = interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{i,j},shape(1),shape(2),'spline');
        end
    end
    
    if isequal(metric,eye(dimension))
        return
    end
    
    %For each shape variable
    for j = 1:dimension

        %Calculate gradient of metric w.r.t. that shape variable using
        %central differencing
        shapeminus = shape;
        shapeminus(j) = shape(j) - shapestep;
        metricminus = zeros(dimension);
        for i = 1:dimension
            for k = 1:dimension
                metricminus(i,k) = interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{i,k},shapeminus(1),shapeminus(2),'spline');
            end
        end

        shapeplus = shape;
        shapeplus(j) = shape(j) + shapestep;
        metricplus = zeros(dimension);
        for i = 1:dimension
            for k = 1:dimension
                metricplus(i,k) = interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{i,k},shapeplus(1),shapeplus(2),'spline');
            end
        end

        thisgrad = (metricplus-metricminus)/(2*shapestep);
        
        %Use formula
        %dMetric/dCoeffs = dMetric/dShape*dShape/dCoeffs
        for i = 1:numel(grad_alpha)
            metricgrad{i} = metricgrad{i}+thisgrad*grad_alpha{i}(j);
        end
        
    end

end