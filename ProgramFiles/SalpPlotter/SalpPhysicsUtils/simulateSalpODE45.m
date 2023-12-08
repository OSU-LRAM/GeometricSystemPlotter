%Evaluates salp physics using ODE45
function [ts,simulationResults] = simulateSalpODE45(salp,forcedef,duration)

    %Set up initial conditions from inputs
    tspan = [0,duration];
    y0 = [salp.geometry.initialPose(:);salp.geometry.initialShape(:)];
    odefun = @(t,y) salpODE(t,y,salp,forcedef);

    %Let the user know how long we expect this to take
    wb = waitbar(0,'Simulating Salp Motion');
    outfun = @(t,y,flag) updateWaitbar(t,y,flag,wb);

    options = odeset('OutputFcn',outfun);

    %Used for updating the waitbar in the middle of ode45
    global totalSimTime; %Total duration of the simulated experiment
    global lastUpdateTime; %Timestamp we last updated the waitbar

    lastUpdateTime = 0;
    totalSimTime = duration;

    %Start timer for waitbar shenanigans
    tic;
    %Run the simulated experiment
    sol = ode45(odefun,tspan,y0,options);

    %Get rid of the waitbar now the experiment is complete
    delete(wb);

    %Resample results based on the framerate we want to make animations in
    FPS = salp.visualization.FrameRate;
    dt = 1/FPS;
    ts = [0:dt:duration];
    simulationResults = deval(sol,ts);

end

%Salp ODE evaluator
function dY = salpODE(t,y,salp,forcedef)

    %[x;y;theta] pose of salp in SE(2)
    pose = y(1:3);
    %value for each shape variable
    shape = y(4:end);

    %Get current thruster force
    f_vec = getForceVectors(salp,forcedef,t);
    %Get velocities that result from that thruster force
    [g_circ,dShape] = salp.ForceConnection(shape,f_vec);

    %Get map from body velocity to world velocity
    R = vec_to_mat_SE2([0,0,pose(3)]);
    %Find world velocity
    dPose = R*g_circ;
    
    %Return state rates of change
    dY = [dPose;dShape];

end

%function to handle waitbar updates while simulation is happening
function haltSim = updateWaitbar(t,y,flag,f)
    
    global totalSimTime;
    global lastUpdateTime;

    %Simulation is currently in the midst of evaluating
    if isempty(flag)

        %What time is sim currently at
        currentTime = t(end);
        %What percent through simulation are we rounded to 1%
        timeFactor = floor(currentTime/totalSimTime * 100)/100;

        %How long has it been since we started simulation
        compTime = toc;
        %If it's been more than a half second since we updated the waitbar
        if (compTime-0.5) > lastUpdateTime

            %Save this time for later updates
            lastUpdateTime = compTime;

            %Estimate how much time is left
            %NOTE: This isn't super accurate, because the simulation is
            %much slower before it hits steady-state.  So this estimate
            %always thinks it's going to take longer than it really will
            ETA = compTime/(currentTime/totalSimTime) - compTime;
            %If time can be measured in seconds do that
            if ETA <= 60
                ETA = floor(ETA);
                waitbar(timeFactor,f,{'Simulating Salp Motion',['Estimated Time Remaining: ',num2str(ETA),' seconds']});
            %Otherwise let the user know it's gonna take fuckin forever
            else
                ETA = floor(ETA/60);
                ETA = ETA + 1;
                waitbar(timeFactor,f,{'Simulating Salp Motion',['Estimated Time Remaining: ',num2str(ETA),' minutes']});
            end
        end

    end

    haltSim = 0;

end


