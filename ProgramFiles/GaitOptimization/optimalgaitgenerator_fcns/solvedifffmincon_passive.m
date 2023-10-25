function [f,g]=solvedifffmincon_passive(y,s,n,dimension,direction,~,~,~)%,lb,ub,writerObj)
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
coeff=[y,y];

if strcmpi(s.passiveObjectiveFunction,'normalized speed')
    coeff(end,:) = [2*pi,2*pi];
end
w1 = coeff(end,2); % Frequency of Fourier transform
% Assign a time period for executing the gait
T = 2*pi/w1;

global gradCoeffList;

if nargout == 1
   
    thisClosestGrad = struct();
    thisClosestDistance = 1e6;
    for i = 1:numel(gradCoeffList)
        coeffs = gradCoeffList{i}.coeffs;
        thisDist = norm(coeff(:,2) - coeffs(:,2));
        if thisDist < thisClosestDistance
            thisClosestDistance = thisDist;
            thisClosestGrad = gradCoeffList{i};
        end
    end
    
    %disp(thisClosestDistance);
    thisDist = coeff(:,2) - thisClosestGrad.coeffs(:,2);
    transferGrads = thisClosestGrad.tf;
    
    newCoeff = thisClosestGrad.coeffs;
    for i = 1:numel(thisDist)
        deltaCoeff = transferGrads{i}*thisDist(i);
        newCoeff = newCoeff + deltaCoeff;
    end
    coeff = newCoeff;
    
else

    %Extract active joint gait as a function of time from fourier coefficients,
    %reframe into functions that the simulator expects
    p_temp = makeGait(y);
    p_active.rc = p_temp.phi_def{1};
    p_active.drc = p_temp.dphi_def{1};
    p_active.ddrc = p_temp.ddphi_def{1};

    %Simulate passive joint behavior, and fit to a fourier series for future
    %estimation of the passive/active transfer function
    [~,~,angles,~] = simulate2DPassiveSwimmer(p_active,T,s.funs,s.physics.k,s.physics.b,0);
    ts = linspace(0,T,numel(angles(1,:)));

    if any(isnan(angles(1,:)))
        angles(1,:) = ts;
    end

    fitOpts = fitoptions('fourier4','StartPoint',y,'Lower',[0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,w1],'Upper',[0,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,w1]);
    passive_fit = fit(ts',angles(1,:)','fourier4',fitOpts);
    nu={'a0';'a1';'b1';'a2';'b2';'a3';'b3';'a4';'b4';'w'};
    for i=1:length(nu)
        coeff(i,1)=passive_fit.(nu{i});
    end

    transferGrads = getTransferGradients(coeff(:,2),coeff(:,1));
    
    gradStruct = struct();
    gradStruct.coeffs = coeff;
    gradStruct.tf = transferGrads;
    gradCoeffList{end+1} = gradStruct;
    
end

global lastCoeffSet lastTransferGrads;
lastCoeffSet = coeff;
lastTransferGrads = transferGrads;

y1 = path_from_fourier(coeff,n,2);
y1 = y1(1:end-1,:); % Remove the last points because path_from_fourier returns self-connected gait



%% Calculating cost and displacement per gait




% Define phi_def = [alpha1, alpha2] as a function of time t such that the
% array returns the shape variables given by the fourier coefficients at a
% time t
p = makeGait(coeff);
        
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

global bestCost bestDisp bestEff bestCoeffs;

% Calculate displacement, cost and efficiency of a gait
% Note that, for inertial cost, cost is returned as the integral of torque
% squared, while for drag-based systems, cost is the path length of the
% gait
[~, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');
lineint=net_disp_opt(direction); % displacement produced in the chosen direction produced on executing the gait measured in the optimal coordinates 

%Group cost types
inertialCosts = {'torque','mechanical power','covariant acceleration','acceleration coord','power quality'};

totalstroke = cost;

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

if nargout > 1

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
            tj = t(j);
            w = coeff(end,i);
            dpdw = -coeff(2,i)*tj*sin(w*tj) + coeff(3,i)*tj*cos(w*tj)...
                -2*coeff(4,i)*tj*sin(2*w*tj) + 2*coeff(5,i)*tj*cos(2*w*tj)...
                -3*coeff(6,i)*tj*sin(3*w*tj) + 3*coeff(7,i)*tj*cos(3*w*tj)...
                -4*coeff(8,i)*tj*sin(4*w*tj) + coeff(9,i)*tj*cos(4*w*tj);
            chy{i}(:,j)=[1;cos(tj*w);sin(tj*w);cos(2*tj*w);sin(2*tj*w);cos(3*tj*w);sin(3*tj*w);cos(4*tj*w);sin(4*tj*w);dpdw];%cos(5*t(j)*coeff(end,i));sin(5*t(j)*coeff(end,i))];%;cos(6*t(j)*coeff(end,i));sin(6*t(j)*coeff(end,i))];%
        end
    end

    %% Jacobianstroke is the gradient of cost. 
    %Contrigrad is the contribution to the gradient due to the metric changing

    switch s.costfunction
        case {'pathlength metric','pathlength coord','pathlength metric2'}
            % Get the gradient of cost based on drag-dominated system
            jacobianstroke = jacobianstrokecalculator(y,n,dimension,metric,metricgrad);
        case 'mechanical power'
            % Get the gradient of cost based on inertia-dominated system
            inertia_cost_grad = inertia_cost_gradient_passive(s,n,coeff,T,p,'discrete');
            if strcmpi(s.passiveObjectiveFunction,'normalized speed')
                totalstroke = totalstroke^(1/3);
                %Rescale gradient to the power-normalized version
                inertia_cost_grad = inertia_cost_grad./(3*totalstroke^2);
            end
        case {'torque','covariant acceleration','acceleration coord','power quality'}
            % Get the gradient of cost based on inertia-dominated system
            inertia_cost_grad = inertia_cost_gradient_passive(s,n,coeff,T,p,'discrete');
            if strcmpi(s.passiveObjectiveFunction,'normalized speed')
                totalstroke = totalstroke^(1/4);
                inertia_cost_grad = inertia_cost_grad./(4*totalstroke^3);
            end
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
        for j=1:1:10
            jacobiandispfourier(j,i)=chy{i}(j,:)*jacobiandisp(:,i);
            if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
                jacobianstrokefourier(j,i)=chy{i}(j,:)*jacobianstroke(:,i);
            end
            jacobianeqifourier(j,i)=chy{i}(j,:)*jacobianeqi(:,i);
        end
    end



    if any(strcmpi(s.costfunction,inertialCosts))
        % Inertia cost gradient is already in terms of the fourier coefficients
        jacobianstrokefourier = inertia_cost_grad;
    end

    costGradHolder = zeros(size(coeff(:,2)));
    distGradHolder = zeros(size(coeff(:,2)));
    equiGradHolder = zeros(size(coeff(:,2)));
    for i = 1:numel(costGradHolder)
        costGradHolder(i) = sum(jacobianstrokefourier.*transferGrads{i},'all');
        distGradHolder(i) = sum(jacobiandispfourier.*transferGrads{i},'all');
        equiGradHolder(i) = sum(jacobianeqifourier.*transferGrads{i},'all');
    end
    jacobianstrokefourier = costGradHolder;
    jacobiandispfourier = distGradHolder;
    jacobianeqifourier = equiGradHolder;

    periodfouriergrad = zeros(size(costGradHolder));
    periodfouriergrad(end) = -2*pi/w1^2;

    switch s.passiveObjectiveFunction
        case 'normalized speed'
            totaljacobianfourier = jacobiandispfourier/totalstroke - lineint*jacobianstrokefourier/totalstroke^2;
            f = -lineint/totalstroke;
        case 'speed'
            totaljacobianfourier = jacobiandispfourier/T - lineint*periodfouriergrad/T^2;
            f = -lineint/T;
        case 'efficiency'
            totaljacobianfourier = jacobiandispfourier/totalstroke - lineint*jacobianstrokefourier/totalstroke^2;
            f = -lineint/totalstroke;
        case 'metabolism'
            metab = s.physics.metabolicRate;
            totaljacobianfourier = jacobiandispfourier/(totalstroke + metab*T)...
                - lineint*(jacobianstrokefourier + metab*periodfouriergrad)/(totalstroke + metab*T)^2;
            f = -lineint/(totalstroke + metab*T);
        otherwise
            error('No passive objective function specified in the system file.  Please choose between "speed", "efficiency", or "metabolism"');
    end

    if ~any(strcmpi(s.costfunction,inertialCosts))
        totaljacobianfourier = totaljacobianfourier+jacobianeqifourier;
    end
    
else

   switch s.passiveObjectiveFunction
       case 'normalized speed'
           if strcmpi(s.costfunction,'mechanical power')
               totalstroke = totalstroke^(1/3);
           elseif strcmpi(s.costfunction,'torque')
               totalstroke = totalstroke^(1/4);
           end
           f = -lineint/totalstroke;
       case 'speed'
           f = -lineint/T;
       case 'efficiency'
           f = -lineint/totalstroke;
       case 'metabolism'
           metab = s.physics.metabolicRate;
           f = -lineint/(totalstroke + metab*T);
       otherwise
           error('No passive objective function specified in the system file.  Please choose between "speed", "efficiency", or "metabolism"');
    end
    
end

if abs(f) > bestEff
    bestEff = abs(f);
    bestDisp = abs(lineint);
    bestCoeffs = coeff(:,2);
    switch s.passiveObjectiveFunction
        case 'normalized speed'
            bestCost = totalstroke;
        case 'speed'
            bestCost = T;
        case 'efficiency'
            bestCost = abs(totalstroke);
        case 'metabolism'
            bestCost = abs(totalstroke + metab*T);
    end
end

printOut = 0;
if printOut
    disp(['Coeffs: ',num2str(coeff(1:5,2)'),' ',num2str(coeff(end,2))]);
    disp(['Distance: ',num2str(abs(lineint))]);
    disp(['Time: ',num2str(T)]);
end

if nargout>1
    
    g = -totaljacobianfourier;
    
    testGrad = 0;
    if testGrad
        oldY = coeff(:,2);
        step = 0.1;
        numericG = zeros(size(g));
        for i = 1:numel(g)
            ytemp = oldY;
            ytemp(i) = oldY(i) + step;
            [newDisp,newT] = quickSim(ytemp,s);
            newF = -abs(newDisp/newT);
            numericG(i) = (newF-f)/step;
        end
        disp(g);
        disp(numericG);
    end
    
end

end