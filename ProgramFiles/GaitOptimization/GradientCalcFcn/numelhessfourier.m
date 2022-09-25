function hess = numelhessfourier(f,s,n,dimension,direction,lb,ub,grad)
dx = 1e-6;
m = size(f,1)-1;
djacobian = cell(m,dimension);
dgrad = cell(m,dimension);
hess = struct();
chy=chy_generator(f,n,dimension);
parfor i = 1:m
    for j = 1:dimension
        df=zeros(m+1,dimension);
        dgrad{i,j} = struct();
        df(i,j) = dx;
        djacobian{i,j} = evaluate_jacobian(f+df,s,n,dimension,direction,lb,ub);
        for k = 1:dimension
            dgrad{i,j}.disp(:,k)=chy{k}*djacobian{i,j}.disp(:,k);
            dgrad{i,j}.stroke(:,k)=chy{k}*djacobian{i,j}.stroke(:,k);
            dgrad{i,j}.repuls(:,k)=chy{k}*djacobian{i,j}.repuls(:,k);
        end
    end
end
hess.disp = cell(dimension,dimension);
hess.stroke = cell(dimension,dimension);
hess.repuls = cell(dimension,dimension);
for i = 1:dimension
    for k = 1:m
        temphessdisp=(dgrad{k,i}.disp - grad.disp).'/dx;
        temphessstroke=(dgrad{k,i}.stroke - grad.stroke).'/dx;
        temphesssrepuls=(dgrad{k,i}.repuls - grad.repuls).'/dx;
        for j = 1:dimension
            hess.disp{i,j}(k,:) = temphessdisp(j,:);
            hess.stroke{i,j}(k,:) = temphessstroke(j,:);
            hess.repuls{i,j}(k,:)= temphesssrepuls(j,:);
        end
    end
end

