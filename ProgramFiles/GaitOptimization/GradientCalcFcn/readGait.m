%Finds value of gait cell function at given time
function state = readGait(gaitFun,t)

    ndim = numel(gaitFun);
    state = zeros(1,ndim);
    
    for i = 1:ndim
        state(i) = gaitFun{i}(t);
    end

end
