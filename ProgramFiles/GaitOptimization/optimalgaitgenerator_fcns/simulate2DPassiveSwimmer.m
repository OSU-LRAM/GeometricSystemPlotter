%Simulates the inertial-passive gait for the swimming system, with the aim
%of finding the limit cycle

%p-position/velocity/acceleration of active joint as function wrt time
%T-period of gait
%funs-metric/coriolis/connection matrices as functions wrt shape
%k-spring constant
%b-damping constant
%animate-boolean, whether or not to make a gif of the resulting motion

%displ-displacement due to limit cycle gait
%cost-torque squared cost to execute limit cycle gait
%angles-joint deflection profile of limit cycle gait
%final_loop-[X motion;Y motion;Theta motion;passive joint;controlled
%joint;cost] over course of limit cycle gait

function [displ,cost,angles,final_loop] = simulate2DPassiveSwimmer(p,T,funs,k,b,animate,loopPassive)

%Store inputs into system structure
sys.k = k;
sys.d = b;
sys.A = funs.A_fun;
sys.dMdr1 = funs.dmdr1_fun;
sys.dMdr2 = funs.dmdr2_fun;
sys.metric = funs.metric_fun;
sys.T = T;
sys.animate_T = 5;

%Initial values for gait
init = [0;0;0;0;0;0];

%Number of cycles to run before evaluating limit cycle
t_lim = 5*T;

%Make gait reader compatible with sysplotter
if ~isfield(p,'rc') && isfield(p,'phi_def')
    p.rc = p.phi_def{end};
    p.drc = p.dphi_def{end};
    p.ddrc = p.ddphi_def{end};
end
sys.p = p;

%Turn animation off by default
if nargin == 5
    animate = 0;
end

%Turn looping of passive joint on by default
if nargin == 6
    loopPassive = 1;
end

sys.loop = 1;

%Set function to evaluate gait dynamics ODE
odefun = @(t,X) getDispAndCost(t,X,sys);
%Evaluate ODE
sol = ode45(odefun,[0,t_lim],init);

%Make even timesteps across limit cycle
nsteps = 100;
dt = T/nsteps;
t = [0:dt:T];

%Initial values of limit cycle loop are ending values of initial loops with
%zeroed displacements and cost
final_loop_start = sol.y(:,end);
final_loop_start(1:3) = [0;0;0];
final_loop_start(end) = 0;

%Evaluate ODE for limit cycle
costDispCalc = @(t,X) getDispAndCost(t,X,sys);
sol_final = ode45(costDispCalc,[0,T],final_loop_start);
sol_final = deval(sol_final,t);

%ODE outputs for outside analysis
final_loop = sol_final;

%Store gait shape positions and velocities
angles = zeros(4,numel(final_loop(1,:)));
angles(1,:) = final_loop(4,:);
angles(2,:) = sys.p.rc(t);
angles(3,:) = final_loop(5,:);
angles(4,:) = sys.p.drc(t);

%Calculate total displacements and cost
gx = final_loop(1,end);
gy = final_loop(2,end);
displ = sqrt(gx^2+gy^2);
cost = final_loop(6,end);

final_loop = addEnergies(t,final_loop,sys);

%If making a gif of output motion, do it
if animate
    t_span = [0,sys.animate_T];
    sol = ode45(odefun,[0,t_span(2)],init);
    t = [0:1/30:t_span(2)];
    sol_a = deval(sol,t);  
    animateSwimmer(sol_a,sys);
end

end

function augmentedStates = addEnergies(ts,Xs,sys)
    
    energies = zeros(3,numel(ts));
    
    for i = 1:numel(ts)
        
        t = ts(i);
        %Get current position/velocity of passive joint
        r1 = Xs(4,i);
        dr1 = Xs(5,i);
        %Get current position/velocity/acceleration of control joint
        r2 = sys.p.rc(t);
        dr2 = sys.p.drc(t);
        ddr2 = sys.p.ddrc(t);
        

        %Wrap joints to pi so they stay in bounds of interpolant functions
        if sys.loop
            r1 = wrapToPi(r1);
        end
        r2 = wrapToPi(r2);
    
        %Interpolate metric from shape
        M = sys.metric(r1,r2);
        
        dr = [dr1;dr2];
        energies(1,i) = .5*dr'*M*dr;
        energies(2,i) = .5*sys.k*r2*r2;
        
        %Get torque on passive joint from passive components
        tau_1 = -sys.k*r1 - sys.d*dr1;
        %Get coriolis forces from shape and shape velocity
        C = getCoriolisForces(r1,dr1,r2,dr2,sys);
        %Calculate acceleration of passive joint using dynamic equation
        ddr1 = (tau_1 - M(1,2)*ddr2 - C(1))/M(1,1);
        %Get required torque at control joint to execute
        tau_2 = M(2,:)*[ddr1;ddr2] + C(2);
        
        energies(3,i) = tau_2;

    end
    
    augmentedStates = [Xs;energies];
    
end
%ODE descriptor for swimmer dynamics with cost included
function dF = getDispAndCost(t,X,sys)    

    dF = zeros(6,1);
    
    %Get current position/velocity of passive joint
    r1 = X(4);
    dr1 = X(5);
    %Get current position/velocity/acceleration of control joint
    r2 = sys.p.rc(t);
    dr2 = sys.p.drc(t);
    ddr2 = sys.p.ddrc(t);
    
    %Get torque on passive joint from passive components
    tau_1 = -sys.k*r1 - sys.d*dr1;

    %Wrap joints to pi so they stay in bounds of interpolant functions
    if sys.loop
        r1 = wrapToPi(r1);
    end
    r2 = wrapToPi(r2);
    
    %Interpolate metric from shape
    M = sys.metric(r1,r2);
    %Get coriolis forces from shape and shape velocity
    C = getCoriolisForces(r1,dr1,r2,dr2,sys);
    
    %Calculate acceleration of passive joint using dynamic equation
    ddr1 = (tau_1 - M(1,2)*ddr2 - C(1))/M(1,1);
    %Get required torque at control joint to execute
    tau_2 = M(2,:)*[ddr1;ddr2] + C(2);
    
    %Interpolate connection
    A = sys.A(r1,r2);
    %Reconstruction equation to get swimmer motion
    g_circ = -A*[dr1;dr2];
    
    %Convert swimmer velocity from local to world coordinates
    theta = X(3);
    g = [cos(theta),-sin(theta),0;...
        sin(theta),cos(theta),0;...
        0,0,1];
    g_dot = g*g_circ;
    
    %Store rate of change of each parameter in ODE
    %World position velocities
    dF(1:3) = g_dot;
    %Passive joint velocity/acceleration
    dF(4:5) = [dr1;ddr1];
    %Cost velocity
    %dF(6) = abs(tau_2*dr2);
    %dF(6) = max(tau_2*dr2,0);
    dF(6) = tau_2*dr2;
    
end

%Calculates coriolis forces from system shape/velocity and derivatives of
%metric
function C = getCoriolisForces(r1,dr1,r2,dr2,sys)

    dr = [dr1;dr2];
    
    C1 = sys.dMdr1(r1,r2)*dr1*dr + sys.dMdr2(r1,r2)*dr2*dr;
    C2 = -.5*[dr'*sys.dMdr1(r1,r2)*dr;dr'*sys.dMdr2(r1,r2)*dr];
    
    C = C1 + C2;
end