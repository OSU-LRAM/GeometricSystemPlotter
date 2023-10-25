function y=passiveoptimalgaitgenerator(s,dimension,npoints,a,lb,ub,stretch,direction,costfunction,handles)
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

% Arbitrarily fit for a starting frequency, time period will be optimized later
startHz = 1.7;
t = linspace(0,1/startHz,size(P1,1));



% The bounds ensure the fourier series terms have the right period
 % Bound only the frequency to be 2*pi, such that period = 1
options = fitoptions('fourier4');
options.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 2*pi*startHz];
options.Upper = -[-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -2*pi*startHz];

fa=fit(t',P1(:,1),'fourier4',options);

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

 lim = 2;
 lb1=[0,0.02,0,0,0,-lim,-lim,0,0,0.5*2*pi];
 ub1=[0,lim,0,0,0,lim,lim,0,0,20];

 y0 = zeros(length(nu),1);
for j=1:length(nu)
    y0(j,1)=fa.(nu{j});
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
s = fitConnectionAndMetric(s);

global bestCost bestDisp bestEff bestCoeffs gradCoeffList;
bestEff = 0;
gradCoeffList = {};

%Suppress warning for annoying thing in jacobianeqicalculator
warning('off','MATLAB:sqrtm:SingularMatrix');

try
 options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
     'Display','off','Algorithm','sqp','CheckGradients',false,'FiniteDifferenceType','central','MaxIter',4000,...
     'MaxFunEvals',20000,'TolCon',10^-4,'StepTolerance',1e-8,...
     'OutputFcn', @(y,optimValues,state) outfun_passive(y,optimValues,state,stretch,s,handles));
catch
   error('This code requires the global optimization toolbox to run') 
end

numRuns = 3;
freqSteps = linspace(1,2,numRuns);

global bestCost bestDisp bestEff bestCoeffs gradCoeffList;

objective_function_gradient = @(y) solvedifffmincon_passive(y,s,npoints,dimension,direction,lb,ub,writerObj);
constraint_function = @(y) nonlcon_passive(y,s,npoints,dimension,lb,ub,direction);

multiStartBestCoeffs = zeros(10,2);
bestScore = 0;
scores = zeros(10,1);
allCoeffs = cell(10,1);
accels = zeros(10,1);

for i = 1:numRuns
    
    disp(num2str(i/numRuns))
    
    w = freqSteps(i);
    y0 = zeros(10,1);
    y0(2) = 0.9*s.physics.maxAcc/(2*pi*w)^2;
    y0(end) = 2*pi*w;
   
    bestEff = 0;
    gradCoeffList = {};

    [yf, ~,~,~]=fmincon(objective_function_gradient,y0,A,b,Aeq,beq,lb1,ub1,constraint_function,options);

    % % Uncomment this if you uncommented the section above so that the video
    % % writer object is closed appropriately.
    % close(writerObj);

    %% Getting point position values from the result of fmincon
    % This section helps us go back to a direct transcription parametrization
    % of the optimal gait from a fourier series parametrization. y is a column vector
    % that contains coordinates of all points forming the optimized gait

    w = yf(end,1);
    T = 2*pi/w;
    coeff=[yf,yf];
    p_temp = makeGait(yf);
    p_active.rc = p_temp.phi_def{1};
    p_active.drc = p_temp.dphi_def{1};
    p_active.ddrc = p_temp.ddphi_def{1};
    [displ,cost,angles,~] = simulate2DPassiveSwimmer(p_active,T,s.funs,s.physics.k,s.physics.b,0);
    ts = linspace(0,T,numel(angles(1,:)));
    passive_fit = fit(ts',angles(1,:)','fourier4','StartPoint',yf);
    nu={'a0';'a1';'b1';'a2';'b2';'a3';'b3';'a4';'b4';'w'};
    for j=1:length(nu)
        coeff(j,1)=passive_fit.(nu{j});
    end

    p = makeGait(coeff);
    [~, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');
    lineint=abs(net_disp_opt(direction)); % displacement produced in the chosen direction produced on executing the gait measured in the optimal coordinates 
    totalstroke = cost;
    
    switch s.passiveObjectiveFunction
        case 'speed'
            f = lineint/T;
        case 'efficiency'
            f = lineint/totalstroke;
        case 'metabolism'
            metab = s.physics.metabolicRate;
            f = lineint/(totalstroke + metab*T);
        otherwise
            error('No passive objective function specified in the system file.  Please choose between "speed", "efficiency", or "metabolism"');
    end
    
    if f > bestScore
        bestScore = f;
        multiStartBestCoeffs = coeff;
    end
    
    scores(i) = f;
    accels(i) = max(abs(p_active.ddrc(ts)));
    allCoeffs{i} = coeff;
    
end

printstuff = 1;
if printstuff
    disp(['Optimal Coeffs: ',num2str(multiStartBestCoeffs(:,2)')]);
    disp(['Optimal Score: ',num2str(bestScore)]);
end

y1 = path_from_fourier(multiStartBestCoeffs,npoints,2);
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