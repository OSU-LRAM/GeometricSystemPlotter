function [metric,metricgrad] = getMetricGrad(s,shape,grad_alpha)
%Returns metric at a shape position, and gradient of metric w.r.t. fourier
%coefficients

%s - structure containing metric function
%shape - array of shape values
%grad_alpha - gradient of shape values w.r.t. fourier coefficients


warning('The function getMetricGrad is deprecated, slow, and should be replaced by interpolation calls to the metric and coriolis fields in the system structure')

    %Number of shape variables
    dimension = numel(shape);
    %Step size for small step approximation on metric gradient
    shapestep = 0.0001;
    
    %Declare empty grid for metric interpolation
    interpmetricgrid=cell(1,dimension);
    for j=1:dimension
        interpmetricgrid{j} = s.grid.metric_eval{j,1};
    end
    
    %Calculate metric at shape position
    shapecell = num2cell(shape);
    
    %Initialize metric gradient to zero matrix for all fourier coefficients
    metricgrad = repmat({zeros(dimension)},size(grad_alpha));
    %If cost function is in coordinate space, metric is always identity
    if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'acceleration coord')
        metric = eye(dimension);
        return
    end
    
    metric = zeros(dimension);
    for i = 1:dimension
        for j = 1:dimension
            metric(i,j) = interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{i,j},shape(1),shape(2),'spline');
        end
    end
    
    if isequal(metric,eye(dimension))
        return
    end
    
    %For each shape variable
    for j = 1:dimension

        %Calculate gradient of metric w.r.t. that shape variable using
        %central differencing
        shapeminus = shape;
        shapeminus(j) = shape(j) - shapestep;
        metricminus = zeros(dimension);
        for i = 1:dimension
            for k = 1:dimension
                metricminus(i,k) = interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{i,k},shapeminus(1),shapeminus(2),'spline');
            end
        end

        shapeplus = shape;
        shapeplus(j) = shape(j) + shapestep;
        metricplus = zeros(dimension);
        for i = 1:dimension
            for k = 1:dimension
                metricplus(i,k) = interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{i,k},shapeplus(1),shapeplus(2),'spline');
            end
        end

        thisgrad = (metricplus-metricminus)/(2*shapestep);
        
        %Use formula
        %dMetric/dCoeffs = dMetric/dShape*dShape/dCoeffs
        for i = 1:numel(grad_alpha)
            metricgrad{i} = metricgrad{i}+thisgrad*grad_alpha{i}(j);
        end
        
    end

end