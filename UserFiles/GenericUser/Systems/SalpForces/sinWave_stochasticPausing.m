function forcedef = sinWave_stochasticPausing(inputmode,nForces,changedDuration)

    forceName = 'Sinusoid Activation - Random Pauses';

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
        
            forcedef.frequency = 0.5;
            forcedef.period = 1/forcedef.frequency;
            forcedef.maxThrust = 1;
            forcedef.nForces = nForces;
        
            %Range of random delay times that will be added between
            %thruster firings
            pauseRange = [0,2];
        
            %Make vector of ignition times that goes up to the force
            %definition duration
            ignitionTimes = cell(1,nForces);
            %For every force
            for i = 1:nForces
                %Start ignition sequence in the past so forces are randomly
                %scattered at the beginning
                thisIgnitionSequence = [-2];
                %Keep adding ignition times until the sequence is long
                %enough
                while(thisIgnitionSequence(end) <= forcedef.duration)
                    thisDelay = pauseRange(1) + rand()*(pauseRange(2)-pauseRange(1));
                    thisIgnitionSequence = [thisIgnitionSequence,thisIgnitionSequence(end)+forcedef.period+thisDelay];
                end
                %Save this ignition sequence to the appropriate thruster
                ignitionTimes{i} = thisIgnitionSequence;
            end
            forcedef.ignitionTimes = ignitionTimes;
        
            %Make cell array of force functions for each thruster
            forceFunctions = cell(forcedef.nForces,1);
            for i = 1:forcedef.nForces
                forceFunctions{i} = @(t) thrustProfile(forcedef,i,t);
            end

            %Call force function by stitching together results of each
            %thruster's individual activation function
            forcedef.amplitudes = @(t) stitchThrustFunctions(forceFunctions,t);
    end
end

%Activation function for this force profile and this thruster
function thrust = thrustProfile(forcedef,thrusterIndex,t)

    %Open ignition sequence for this thruster
    ignitionTimes = forcedef.ignitionTimes{thrusterIndex};
    %Find the ignition index that most recently occurred
    ignitionIndex = find(ignitionTimes>t,1)-1;
    %Get that ignition time
    ignitionTime = ignitionTimes(ignitionIndex);
    %Time since ignition
    dt = t-ignitionTime;

    %Max thrust for scaling
    maxThrust = forcedef.maxThrust;

    %If in the pause after firing set thrust to zero
    if dt > forcedef.period
        thrust = 0;
    %Otherwise set thrust using sinusoid
    else
        thrust = maxThrust/2-maxThrust/2*cos(2*pi*forcedef.frequency*dt);
    end

end