function s = ensure_connection_and_metric(s)


    % Get the number of shape variables
	n_shape = nargin(s.A_num);
	shape_test_list = num2cell(ones(n_shape,1));

    % If there is no connection denominator, make it all ones
    if ~isfield(s,'A_den')
        
		shape_test_list = num2cell(ones(n_shape,1));
		s.A_den = @(a,varargin) repmat(ones(size(s.A_num(shape_test_list{:}))),size(a));
        
        % Indicate that this can be processed as a block-formatted array
        s.function_type.A_den = 'block';
        
    end
    

	%Ensure presence of a metric and a metric denominator
    if ~isfield(s,'metric')
        
		s.metric = @(a,varargin) repmat(eye(size(shape_test_list,1)),size(a));
        
        % Indicate that this should be processed as a woven-formatted array
        s.function_type.metric = 'woven';
        
    end

    if ~isfield(s,'metric_den')
        
		s.metric_den = @(a,varargin) repmat(ones(size(shape_test_list,1)),size(a));
        
        % Indicate that this can be processed as a block-formatted array
        s.function_type.metric_den = 'block';
        
    end
    

    
end