function [jacobianstroke,jacobianeqi] = pathlengthcost(y,M,dM,P)

n = P.N;
dimension = P.shv;
jacobianstroke = zeros(n,dimension);
contrigrad=zeros(n,dimension);

% Calculate the length of segments
for i=1:1:n-1
    L(i) = sqrt((y(i+1,:)-y(i,:))*((M{i}+M{i+1})/2)*(y(i+1,:)-y(i,:))');
end
L(n)=sqrt((y(1,:)-y(n,:))*((M{n}+M{1})/2)*(y(1,:)-y(n,:))');


%% Jacobiandisp-jacobian for displacement produced by gait
for i=1:n-1
    delp{i}=y(i+1,:)-y(i,:);
end
delp{n}=y(1,:)-y(n,:);

for i=2:n-1
    for j=1:dimension
        contrigrad(i,j)=0.5*delp{i}*dM{i,j}*delp{i}'/(2*L(i))+0.5*delp{i-1}*dM{i,j}*delp{i-1}'/(2*L(i-1));
    end
    jacobianstroke(i,:)=(-(((M{i}+M{i+1})/2)*delp{i}')'-(delp{i}*((M{i}+M{i+1})/2)))/(2*L(i))+((((M{i-1}+M{i})/2)*delp{i-1}')'+(delp{i-1}*((M{i}+M{i-1})/2)))/(2*L(i-1))+contrigrad(i,:);%+0.5*[delp{i}*metricgradx{i}*delp{i}',delp{i}*metricgrady{i}*delp{i}',delp{i}*metricgradz{i}*delp{i}']/(2*l(i))+0.5*[delp{i-1}*metricgradx{i}*delp{i-1}',delp{i-1}*metricgrady{i}*delp{i-1}',delp{i-1}*metricgradz{i}*delp{i-1}']/(2*l(i-1));
end

for j=1:dimension
    contrigrad(1,j)=0.5*delp{1}*dM{1,j}*delp{1}'/(2*L(1))+0.5*delp{n}*dM{1,j}*delp{n}'/(2*L(n));
end
jacobianstroke(1,:)=(-(((M{1}+M{2})/2)*delp{1}')'-(delp{1}*((M{1}+M{2})/2)))/(2*L(1))+((((M{n}+M{1})/2)*delp{n}')'+(delp{n}*((M{n}+M{1})/2)))/(2*L(n))+contrigrad(1,:);%+0.5*[delp{1}*metricgradx{1}*delp{1}',delp{1}*metricgrady{1}*delp{1}',delp{1}*metricgradz{1}*delp{1}']/(2*l(1))+0.5*[delp{n}*metricgradx{1}*delp{n}',delp{n}*metricgrady{1}*delp{n}',delp{n}*metricgradz{1}*delp{n}']/(2*l(n));   

for j=1:dimension
    contrigrad(n,j)=0.5*delp{n}*dM{n,j}*delp{n}'/(2*L(n))+0.5*delp{n-1}*dM{n,j}*delp{n-1}'/(2*L(n-1));
end
jacobianstroke(n,:)=(-(((M{n}+M{1})/2)*delp{n}')'-(delp{n}*((M{n}+M{1})/2)))/(2*L(n))+((((M{n}+M{n-1})/2)*delp{n-1}')'+(delp{n-1}*((M{n}+M{n-1})/2)))/(2*L(n-1))+contrigrad(n,:);%+0.5*[delp{n}*metricgradx{n}*delp{n}',delp{n}*metricgrady{n}*delp{n}',delp{n}*metricgradz{n}*delp{n}']/(2*l(n))+0.5*[delp{n-1}*metricgradx{n}*delp{n-1}',delp{n-1}*metricgrady{n}*delp{n-1}',delp{n-1}*metricgradz{n}*delp{n-1}']/(2*l(n-1));   


%% Equi-distance jacobian   
for i=2:n-1
    len=sqrt((y(i+1,:)-y(i-1,:))*((M{i-1}+M{i+1})/2)*(y(i+1,:)-y(i-1,:))');
    midpoint=y(i-1,:)+((y(i+1,:)-y(i-1,:))*sqrtm((M{i-1}+M{i+1})/2))/2;
    betacos=(y(i+1,:)-y(i-1,:))*sqrtm((M{i-1}+M{i+1})/2)*((y(i,:)-y(i-1,:))*sqrtm((M{i-1}+M{i})/2))'/(L(i-1)*len);
    xhat=y(i-1,:)+(y(i+1,:)-y(i-1,:))*sqrtm((M{i-1}+M{i+1})/2)*L(i-1)*betacos/len;
    jacobianeqi(i,:)=midpoint-xhat;
end

    len=sqrt((y(2,:)-y(n,:))*((M{2}+M{n})/2)*(y(2,:)-y(n,:))');
    midpoint=y(n,:)+((y(2,:)-y(n,:))*sqrtm((M{n}+M{2})/2))/2;
    betacos=(y(2,:)-y(n,:))*sqrtm((M{n}+M{2})/2)*((y(1,:)-y(n,:))*sqrtm((M{n}+M{1})/2))'/(L(n)*len);
    xhat=y(n,:)+(y(2,:)-y(n,:))*sqrtm((M{n}+M{2})/2)*L(n)*betacos/len;
    jacobianeqi(1,:)=(midpoint-xhat);

    len=sqrt((y(1,:)-y(n-1,:))*((M{1}+M{n-1})/2)*(y(1,:)-y(n-1,:))');
    midpoint=y(n-1,:)+((y(1,:)-y(n-1,:))*sqrtm((M{1}+M{n-1})/2))/2;
    betacos=(y(1,:)-y(n-1,:))*sqrtm((M{n-1}+M{1})/2)*((y(n,:)-y(n-1,:))*sqrtm((M{n-1}+M{n})/2))'/(L(n-1)*len);
    xhat=y(n-1,:)+(y(1,:)-y(n-1,:))*sqrtm((M{n-1}+M{1})/2)*L(n-1)*betacos/len;
    jacobianeqi(n,:)=(midpoint-xhat);

