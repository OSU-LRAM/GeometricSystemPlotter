function forcedef = simpleSinWave(inputmode,nForces,changedDuration)

    forceName = 'Sinusoidal Activation - No Pausing';

    switch inputmode
        case 'getName'
            forcedef = forceName;
        case 'getForceDefinition'
            forcedef.name = forceName;
            %Default force duration in seconds
            forcedef.duration = 10;
            if changedDuration
                forcedef.duration = changedDuration;
            end
        
            forcedef.phaseOffset = -pi/4;
            forcedef.frequency = 0.5;
            forcedef.maxThrust = 1;
        
            forcedef.nForces = nForces;
        
            forceFunctions = cell(forcedef.nForces,1);
            for i = 1:forcedef.nForces
                forceFunctions{i} = @(t) thrustProfile(forcedef,i,t);
            end
        
            forcedef.amplitudes = @(t) stitchThrustFunctions(forceFunctions,t);
    end
end

%Thrust profile for each thruster for this force definition
function thrust = thrustProfile(forcedef,thrusterIndex,t)

    maxThrust = forcedef.maxThrust;

    thrust = maxThrust/2-maxThrust/2*cos(2*pi*forcedef.frequency*t + (thrusterIndex-1)*forcedef.phaseOffset);

end