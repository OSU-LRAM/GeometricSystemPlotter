%Convert cell array of individual thrust functions into vector of thrust
%values
function allThrusts = stitchThrustFunctions(individualForceFunctions,t)

    nThrusts = numel(individualForceFunctions);
    allThrusts = zeros(nThrusts,1);
    for i = 1:nThrusts
        allThrusts(i) = individualForceFunctions{i}(t);
    end

end