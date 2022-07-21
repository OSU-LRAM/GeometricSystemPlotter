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

% set the number of the fourier coefficient.
nfourier = 4;
strfit = strcat('fourier',num2str(nfourier));

% Time period of gait is 1 second in order to handle calculations performed
% for inertial gait optimization
t = linspace(0,1,size(P1,1));

fa=cell(dimension);
% The bounds ensure the fourier series terms have the right period
% Bound only the frequency to be 2*pi, such that period = 1
options = fitoptions(strfit);
options.Lower = [-Inf(1,2*nfourier + 1) 2*pi];
options.Upper = -[-Inf(1,2*nfourier + 1) -2*pi];

for i=1:1:dimension
    fa{i}=fit(t',P1(:,i),strfit,options);
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

nfparam = (nfourier+1)*2;
nu = cell(nfparam,1);

% generate the coefficient string, based on the parameter number.
% nu={'a0';'a1';'b1';'a2';'b2';'a3';'b3';'a4';'b4';'w'};
for i = 1:nfparam
    if i == 1
        nu{i} = 'a0';
    elseif i == nfparam
        nu{i} = 'w';
    elseif mod(i,2) == 0
        nu{i} = strcat('a',num2str(floor(i/2)));
    elseif mod(i,2) == 1
        nu{i} = strcat('b',num2str(floor(i/2)));
    end
end
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
global bestCost bestDisp bestEff;
bestEff = 0;

%Suppress warning for annoying thing in jacobianeqicalculator
warning('off','MATLAB:sqrtm:SingularMatrix');

try
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
        'Display','iter','Algorithm','sqp','CheckGradients',false, ...
        'FiniteDifferenceType','central','MaxIter',4000, ...
        'MaxFunEvals',20000,'TolCon',10^-3,'StepTolerance',1e-6, ...
        'OutputFcn', @(y,optimValues,state) outfun(y,optimValues,state,stretch,s,handles));
catch
    error('This code requires the global optimization toolbox to run')
end


objective_function_gradient = @(y) solvedifffmincon(y,s,npoints,dimension,direction,lb,ub,writerObj);
constraint_function = @(y) nonlcon(y,s,npoints,dimension,lb,ub,direction);

[yf, ~,~,~]=fmincon(objective_function_gradient,y0,A,b,Aeq,beq,lb1,ub1,constraint_function,options);

% % Uncomment this if you uncommented the section above so that the video
% % writer object is closed appropriately.
% close(writerObj);

printstuff = 1;
if printstuff
    disp(['Optimal Efficiency: ',num2str(bestEff)]);
    disp(['Optimal Displacement: ',num2str(bestDisp)]);
    disp(['Optimal Cost: ',num2str(bestCost)]);
end

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

















