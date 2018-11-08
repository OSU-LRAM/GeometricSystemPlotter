function y=inertial_optimalgaitgenerator(s,dimension,npoints,a1,a2,lb,ub,n_interp)

addpath('Tools')

P.N = npoints;
P.shv = dimension;
P.T = 10;   % Total time
P.dt = P.T/P.N;
%----I believe Hossein said all these P variables were not needed----%
% % Select the method you want to obtain the chrstoffel symbols: "Analytical" or "Numerical"
% P.method = 'Numerical';
% % Select type of cost you want to solve "Delta-v" or "PathLength"
% P.costType = 'Delta-v';
% % define the initial velocity "noncte" or "cte"
% P.initial_velocity = 'noncte';

%--------Start load--------%
% % load('sysf_honey_swimmer_4link_continuous_calc.mat');
% % load('sysf_three_link_lowRe_calc.mat');
% % load('sysf_floating_snake_calc5')
% 
% % load('sysf_floating_snake_calc')
% 
% % load('drag_max_gait')
% % load('sysf_serpenoid_lowRe_calc');
% % sysf_serpenoid_lowRe_calc
% 
% % load('sysf_water_swimmer_calc3')
% 
% load('sysf_threelink_highRe_coupled_calc')
% 
% % load('sysf_serpenoid_highRe_calc')
%--------End load--------%

%coordinates of the points

% for i=1:1:P.N
%     P1(i,1)= 0.9*cos((i-1)*2*pi/P.N);
%     P1(i,2)= 0.9*sin((i-1)*2*pi/P.N);
% end
% P1 = flip(P1);

%-----Substituting the format from optimalgaitgenerator.m-----%
P1(:,1)=a1(1,2:P.N-1)';
P1(:,2)=a2(1,2:P.N-1)';
t1 = linspace(0,1,P.N-2);
t2 = linspace(0,1,n_interp);
P1_new = interp1(t1,P1(:,1),t2,'linear');
P2_new = interp1(t1,P1(:,2),t2,'linear');

% Setting up fmincon

P0=[P1(:,1);P1(:,2)];
P0_new = [P1_new';P2_new'];
P.N = n_interp;
%-----Commenting out because it says it's unused-----%
% time_limit = 1400;
% t0 = [0 time_limit];

% Extract the metric
% M = celltensorconvert(s.metricfield.eval.content.metric);

%-----Commenting this out since it says it's unused-----%
% M = celltensorconvert(s.metricfield.metric_eval.content.metric);
% M = celltensorconvert(M);
% 
% h = metricellipsefield(s.grid.eval{1},s.grid.eval{2},M,'tissot',{});


%% ODE Solver for test
% options = odeset('RelTol',1e-6,'AbsTol',1e-8,'Stats','on','OutputFcn',@odeplot,'MaxStep','1');
% 
% sol = ode45(@(t,Y) optimizer_calc(t,Y,s,P),t0,P0,options);
%%

A=[];
b=[];
Aeq=[];
beq=[];
%-----lb and ub are passed in, commenting out-----%
% lb=-1.7*ones(size(P0));
% ub= 1.7*ones(size(P0));
nonlcon=[];

options = optimoptions('fmincon','Display','iter','Algorithm','active-set','GradObj','on','TolX',1*10^-5,'TolFun',1*10^-5,'MaxIter',2000,'MaxFunEvals',20000);
tic;
[y fval exitflag output] = fmincon(@(Y) inertial_optimizer_calc(Y,s,P),P0_new,A,b,Aeq,beq,lb,ub,nonlcon,options);
toc
% Interpolate back to the expected number of points
t1 = linspace(0,1,npoints);
t2 = linspace(0,1,n_interp);
y1 = interp1(t2,y(1:n_interp),t1);
y2 = interp1(t2,y(n_interp+1:end),t1);
y = [y1 y2]';
%-----Commenting out the plotting stuff-----%
% figure(15)
% hold 'on'
% plot(y(1:n_interp),y(n_interp+1:2*n_interp),'g--');
% axis equal
% 
% pause(0.1)

end





