function metricgrad = getMetricGrad(s,shape,dM,grad_alpha)
%Returns metric at a shape position, and gradient of metric w.r.t. fourier
%coefficients

%s - structure containing metric function
%shape - array of shape values
%dM - derivative of matrix with respect to shape variables
%grad_alpha - gradient of shape values w.r.t. fourier coefficients

    n_dim = numel(s.grid.eval);
    actual_size = min(size(shape,2),n_dim);
    
    metricgrad = repmat({zeros(actual_size)},size(grad_alpha));
    for j = 1:actual_size
        
        for k = 1:numel(grad_alpha)
            metricgrad{k} = metricgrad{k}+ dM{j}*grad_alpha{k}(j);
        end
    end

end