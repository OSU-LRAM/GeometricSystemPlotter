function [coeff,efficiency] = generatePassiveGait_byFreqAndConstraint(s,freq,dimension,npoints,a,lb,ub,stretch,direction,costfunction,y0)
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
startHz = freq;
t = linspace(0,1/startHz,size(P1,1));

% y0 = [0,startAmp,0,0,0,0,0,0,0,2*pi*startHz]';


% The bounds ensure the fourier series terms have the right period
 % Bound only the frequency to be 2*pi, such that period = 1
options = fitoptions('fourier4');
options.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 2*pi*startHz];
options.Upper = -[-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -2*pi*startHz];

%fa=fit(t',P1(:,1),'fourier4',options);

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

if freq

%     lb1=[0,0.01,0,-lim/2,-lim/2,-lim/3,-lim/3,-lim/4,-lim/4,2*pi*startHz];
%     ub1=[0,lim,0,lim/2,lim/2,lim/3,lim/3,lim/4,lim/4,2*pi*startHz];
    
    lb1=[0,0.01,0,0,0,-lim/3,-lim/3,0,0,2*pi*startHz];
    ub1=[0,lim,0,0,0,lim/3,lim/3,0,0,2*pi*startHz];
    
else
    
    lb1=[0,0.01,0,-lim/2,-lim/2,-lim/3,-lim/3,-lim/4,-lim/4,pi];
    ub1=[0,lim,0,lim/2,lim/2,lim/3,lim/3,lim/4,lim/4,4*2*pi];
    
end

%  y0 = zeros(length(nu),1);
% for j=1:length(nu)
%     y0(j,1)=fa.(nu{j});
% end

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
% s = fitConnectionAndMetric(s);

global bestCost bestDisp bestEff bestCoeffs;
bestEff = 0;

%Suppress warning for annoying thing in jacobianeqicalculator
warning('off','MATLAB:sqrtm:SingularMatrix');

try
 options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
     'Display','iter','Algorithm','sqp','CheckGradients',false,'FiniteDifferenceType','central','MaxIter',4000,...
     'MaxFunEvals',20000,'TolCon',10^-5,'StepTolerance',1e-9);
 %,'Algorithm','sqp'
catch
   error('This code requires the global optimization toolbox to run') 
end


objective_function_gradient = @(y) solvedifffmincon_passive(y,s,npoints,dimension,direction,lb,ub,writerObj);
constraint_function = @(y) nonlcon_passive(y,s,npoints,dimension,lb,ub,direction);

[yf,fval]=fmincon(objective_function_gradient,y0,A,b,Aeq,beq,lb1,ub1,constraint_function,options);

% % Uncomment this if you uncommented the section above so that the video
% % writer object is closed appropriately.
% close(writerObj);

printstuff = 0;
if printstuff
    disp(['Optimal Efficiency: ',num2str(bestEff)]);
    disp(['Optimal Displacement: ',num2str(bestDisp)]);
    disp(['Optimal Cost: ',num2str(bestCost)]);
    disp(['Optimal Coeffs: ',num2str(bestCoeffs(:)')]);
end

%% Getting point position values from the result of fmincon
% This section helps us go back to a direct transcription parametrization
% of the optimal gait from a fourier series parametrization. y is a column vector
% that contains coordinates of all points forming the optimized gait

w = yf(end,2);
T = 2*pi/w;
coeff=[yf,yf];
p_temp = makeGait(yf);
p_active.rc = p_temp.phi_def{1};
p_active.drc = p_temp.dphi_def{1};
p_active.ddrc = p_temp.ddphi_def{1};
[displ,cost,angles,~] = simulate2DPassiveSwimmer(p_active,T,s.funs,s.physics.k,s.physics.b,0);
switch s.passiveObjectiveFunction
    case 'normalized speed'
        T = cost^(1/3);
        efficiency = abs(displ/T);
    case 'speed'
        efficiency = abs(displ/T);
end
%disp(['Final Speed: ',num2str(abs(displ/T))]);
fitOpts = fitoptions('fourier4','StartPoint',yf,'Lower',[0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,w],'Upper',[0,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,w]);
ts = linspace(0,T,numel(angles(1,:)));
passive_fit = fit(ts',angles(1,:)','fourier4',fitOpts);
nu={'a0';'a1';'b1';'a2';'b2';'a3';'b3';'a4';'b4';'w'};
for i=1:length(nu)
    coeff(i,1)=passive_fit.(nu{i});
end

y1 = path_from_fourier(coeff,npoints,dimension);
% path_from_fourier returns a self-connected gait, so remove the last point
% to give what optimalgaitgenerator expects to return
y1 = y1(1:end-1,:);
y=y1(:);

%disp('Final Params:')
%disp(yf);

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