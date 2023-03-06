function status = odeoutfun(t,y,flag,s,direction,dimension,nfparam,stretch,handles)
    status = 0;
    
    if strcmpi(flag,'init')
        disp('Step-optimizer starts');
    elseif strcmpi(flag,'done')
        disp('Step-optimizer ends');
    else
        global previousDisp currentDisp currentCost Eff rotOptMod;
        if direction == 4
            dir = 1;
        else
            dir = direction;
        end
        disp(['Current Displacement: ',num2str(currentDisp(dir))]);
        ys=reshape(y(:,end),[nfparam dimension]);
        optimValues.iteration = length(t);
        outfun(ys,optimValues,0,stretch,s,handles);
        
        if(direction == 4)
            % Eff = [Eff; currentDisp(3)/currentCost currentDisp(1)/currentCost];
        else
            Eff = [Eff; currentDisp(direction) currentDisp(direction)/currentCost];
        end
        
        if length(t) > 1
            dt = t(end)-t(end-1);
        else
            dt = t;
        end
        
        %         if abs(currentDisp(dir) - previousDisp) < 1e-5
        if currentDisp(dir) < 1e-2
            previousDisp = 0;
            rotOptMod = 1;
            status=1;
        else
            previousDisp = currentDisp(dir);
        end
    end
end

