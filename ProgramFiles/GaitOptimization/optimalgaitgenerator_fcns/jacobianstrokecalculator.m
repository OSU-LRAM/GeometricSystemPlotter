function jacobianstroke = jacobianstrokecalculator(y,n,dimension,metric,metricgrad)
    % Calculates the gradient of cost for drag dominated systems.
    % Inputs:
    %   y: matrix containing points that transcribe the gait
    %   n: number of points in gait transcription
    %   dimension: number of shape variables
    %   metric: Riemannian metric
    %   metricgrad: Gradient of Riemannian metric
    
    %l is the vector containing metric weighted distances between neighbouring
    %points
    l = zeros(1,n);
    for i=1:(numel(l)-1)
        l(i)=sqrt((y(i+1,:)-y(i,:))*((metric{i}+metric{i+1})/2)*(y(i+1,:)-y(i,:))');
    end
    l(end)=sqrt((y(1,:)-y(n,:))*((metric{n}+metric{1})/2)*(y(1,:)-y(n,:))');
    
    delp = cell(1,n);
    for i=1:(numel(delp)-1)
        delp{i}=y(i+1,:)-y(i,:); % delp{i} is the vector joining the (i+1)th point to the ith point 
    end
    delp{end}=y(1,:)-y(n,:);

    jacobianstroke = zeros(n,dimension);
    contrigrad=zeros(n,dimension);
    for i=2:n-1
        for j=1:dimension
            %Contrigrad is the contribution to the gradient due to the metric changing
            contrigrad(i,j)=0.5*delp{i}*metricgrad{i,j}*delp{i}'/(2*l(i))+0.5*delp{i-1}*metricgrad{i,j}*delp{i-1}'/(2*l(i-1)); 
        end
        % Total gradient is the result of distance changing due to movement of point and the metric changing due to movement of the point
        jacobianstroke(i,:)=(-(((metric{i}+metric{i+1})/2)*delp{i}')'-(delp{i}*((metric{i}+metric{i+1})/2)))/(2*l(i))+...
            +((((metric{i-1}+metric{i})/2)*delp{i-1}')'+(delp{i-1}*((metric{i}+metric{i-1})/2)))/(2*l(i-1))+contrigrad(i,:); 
    end

    % Calculation for the 1st point and last point have to be done outside the
    % loop as the (i+1)th point for the last point is the first point and
    % (i-1)th point for the first point is the last point
    for j=1:dimension
        contrigrad(1,j)=0.5*delp{1}*metricgrad{1,j}*delp{1}'/(2*l(1))+0.5*delp{n}*metricgrad{1,j}*delp{n}'/(2*l(n));
    end
    jacobianstroke(1,:)=(-(((metric{1}+metric{2})/2)*delp{1}')'-(delp{1}*((metric{1}+metric{2})/2)))/(2*l(1))+...
        +((((metric{n}+metric{1})/2)*delp{n}')'+(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+contrigrad(1,:);

    for j=1:dimension
        contrigrad(n,j)=0.5*delp{n}*metricgrad{n,j}*delp{n}'/(2*l(n))+0.5*delp{n-1}*metricgrad{n,j}*delp{n-1}'/(2*l(n-1));
    end
    jacobianstroke(n,:)=(-(((metric{n}+metric{1})/2)*delp{n}')'-(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+...
        +((((metric{n}+metric{n-1})/2)*delp{n-1}')'+(delp{n-1}*((metric{n}+metric{n-1})/2)))/(2*l(n-1))+contrigrad(n,:);
end