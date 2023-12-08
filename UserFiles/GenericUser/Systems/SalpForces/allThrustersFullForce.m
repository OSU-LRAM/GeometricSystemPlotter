function forcedef = allThrustersFullForce(inputmode,nForces,changedDuration)

    %I'm a silly boy
    forceName = 'FULL BLAST BABY';

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
        
            forcedef.maxThrust = 1;
            forcedef.nForces = nForces;
        
            forceFunctions = cell(forcedef.nForces,1);
            for i = 1:forcedef.nForces
                forceFunctions{i} = @(t) thrustProfile(forcedef,i,t);
            end
        
            forcedef.amplitudes = @(t) stitchThrustFunctions(forceFunctions,t);
    end
end

%FULL STEAM AHEAD
function thrust = thrustProfile(forcedef,thrusterIndex,t)

    thrust = 1;

end