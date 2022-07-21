function jacobianeqi = jacobianeqicalculator(y,n,dimension,metric)
    % Calculates the gradient of the force driving the points along the
    % gait to be equally spaced.
    % Inputs:
    %   y: matrix containing points that transcribe the gait
    %   n: number of points in gait transcription
    %   dimension: number of shape variables
    %   metric: Riemannian metric
    
    jacobianeqi = zeros(n,dimension);
    
    %l is the vector containing metric weighted distances between neighbouring
    %points
    l = zeros(1,n);
    for i=1:(numel(l)-1)
        l(i)=sqrt((y(i+1,:)-y(i,:))*((metric{i}+metric{i+1})/2)*(y(i+1,:)-y(i,:))');
    end
    l(end)=sqrt((y(1,:)-y(n,:))*((metric{n}+metric{1})/2)*(y(1,:)-y(n,:))');
    
    for i=2:n-1
        len=sqrt((y(i+1,:)-y(i-1,:))*((metric{i-1}+metric{i+1})/2)*(y(i+1,:)-y(i-1,:))'); % metric weighted length between point (i-1) and (i+1)
        midpoint=y(i-1,:)+((y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2))/2; % location of midpoint of the line segment joining point (i-1) and (i+1)
        betacos=(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*((y(i,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i})/2))'/(l(i-1)*len);
        xhat=y(i-1,:)+(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*l(i-1)*betacos/len; %projection of ith point onto the line joining the (i-1)th and (i+1)th points
        jacobianeqi(i,:)=midpoint-xhat; % gradient of the ith point is equal to the difference between the midpoint and the projection of ith point
    end

    len=sqrt((y(2,:)-y(n,:))*((metric{2}+metric{n})/2)*(y(2,:)-y(n,:))');
    midpoint=y(n,:)+((y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2))/2;
    betacos=(y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2)*((y(1,:)-y(n,:))*sqrtm((metric{n}+metric{1})/2))'/(l(n)*len);
    xhat=y(n,:)+(y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2)*l(n)*betacos/len;
    jacobianeqi(1,:)=midpoint-xhat;

    len=sqrt((y(1,:)-y(n-1,:))*((metric{1}+metric{n-1})/2)*(y(1,:)-y(n-1,:))');
    midpoint=y(n-1,:)+((y(1,:)-y(n-1,:))*sqrtm((metric{1}+metric{n-1})/2))/2;
    betacos=(y(1,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{1})/2)*((y(n,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{n})/2))'/(l(n-1)*len);
    xhat=y(n-1,:)+(y(1,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{1})/2)*l(n-1)*betacos/len;
    jacobianeqi(n,:)=midpoint-xhat;
end