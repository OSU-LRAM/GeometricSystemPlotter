%DEPRECATED
%This was an attempt to use RK4 euler's method to integrate salp physics
%Turns out that it's both less accurate than ODE45 and requires setting an
%appropriate time-step to get a stable solution that doesn't take forever
%So this probably won't be used anymore
function [ts,simulationResults] = simulateSalpExponentialEuler(salp,forcedef,duration)

    dt = 1/500;
    ts = [0:dt:duration];
    
    simResults = zeros(3+salp.geometry.nShapes,numel(ts));
    simResults(:,1) = [salp.geometry.initialPose(:);salp.geometry.initialShape(:)];

    wb = waitbar(0,'Simulating Salp Motion');

    completeFactor = 0;
    for i = 1:numel(ts)-1

        y = simResults(:,i);
        t = ts(i);
        simResults(:,i+1) = exponentialEulerTimeStepRK4(t,y,salp,forcedef,dt);

        timeFactor = floor(t/duration*100)/100;
        if timeFactor > completeFactor
            completeFactor = timeFactor;
            waitbar(timeFactor,wb);
        end

    end

    delete(wb);

    FPS = salp.visualization.FrameRate;
    newdt = 1/FPS;
    newts = [0:newdt:duration];

    simulationResults = interp1(ts,simResults',newts)';
    ts = newts;

end

function nextState = exponentialEulerTimeStepRK4(t,y,salp,forcedef,dt)

    k1 = salpStateDerivative(t,y,salp,forcedef);
    y1 = exponentialIntegrateState(y,k1,dt/2);
    k2 = salpStateDerivative(t+dt/2,y1,salp,forcedef);
    y2 = exponentialIntegrateState(y,k2,dt/2);
    k3 = salpStateDerivative(t+dt/2,y2,salp,forcedef);
    y3 = exponentialIntegrateState(y,k3,dt);
    k4 = salpStateDerivative(t+dt,y3,salp,forcedef);

    dy = (k1 + 2*k2 + 2*k3 + k4)/6;
    nextState = exponentialIntegrateState(y,dy,dt);

end

function dy = salpStateDerivative(t,y,salp,forcedef)

    pose = y(1:3);
    shape = y(4:end);

    f_vec = getForceVectors(salp,forcedef,t);
    [g_circ,dShape] = salp.ForceConnection(shape,f_vec);

    R = vec_to_mat_SE2([0,0,pose(3)]);
    dPose = R*g_circ;

    dy = [dPose(:);dShape(:)];

end

function newState = exponentialIntegrateState(y,dy,dt)

    pose = y(1:3);
    shape = y(4:end);

    dpose = dy(1:3);
    dshape = dy(4:end);

    poseMatrix = vec_to_mat_SE2(pose);
    dposeMatrix = vec_to_mat_SE2_lie(dpose);
    newPoseMatrix = poseMatrix*expm(dt*dposeMatrix);
    newPose = mat_to_vec_SE2(newPoseMatrix);

    newShape = shape + dt*dshape;

    newState = [newPose(:);newShape(:)];

end