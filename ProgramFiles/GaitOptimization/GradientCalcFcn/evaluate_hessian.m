function hessian = evaluate_hessian(y,s,n,dimension,direction,invert)

if invert == 1
    y=flip(y);
end

%% Hessian Calculator
hessian = struct();
hessian.disp=hessiandispcalculator3(y,s,n,dimension,direction);
hessian.stroke=numelhessstroke(y,s,n,dimension);

if invert == 1
    for i = 1:dimension
        for j = 1:dimension
            hessian.disp{i,j}=flip(flip(hessian.disp{i,j},1),2);
            hessian.stroke{i,j}=flip(flip(hessian.stroke{i,j},1),2);
        end
    end
end
end