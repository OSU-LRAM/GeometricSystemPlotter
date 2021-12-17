% ensure_connection_and_metric.m
% sanitizes local connection, metric, and group definition for a system

function s = ensure_connection_and_metric(s)
    %% connection checks
    % If the connection was specified as s.A, copy this into s.A_num
    if ~isfield(s,'A_num')
        if isfield(s,'A')
            s.A_num = s.A;
        else
            error('No local connection was specified as s.A or s.A_num in the sysf_file')
        end
    end
    % Get the number of shape variables
	n_shape = nargin(s.A_num);
    
    % Choose a random point within the plot range (to probabalistically
    % avoid testing at a singularity)
    shape_test_point = ones(n_shape,1);
    for idx = 1:n_shape
       shape_test_point(idx) = ...
           (rand(1)*(s.grid_range(2*idx)-s.grid_range((2*idx)-1)))-s.grid_range((2*idx)-1);
    end
	shape_test_list = num2cell(shape_test_point);
    
    % Ensure existence of configuration space
    if ~isfield(s, 'conf_space')
        s.conf_space = LieGroups.SE2;
        warning('No configuration space specified; assuming SE(2). To avoid this warning, add "s.conf_space = LieGroups.SE2" to the sysf.');
    end
    % Enforce configuration space dimensions
    A_dim = size(s.A_num(shape_test_list{:}),1);
    if A_dim ~= s.conf_space.n_dim
        error(['Supplied local connection is not of same dimensionality (' num2str(A_dim) ') as configuration space (' num2str(s.conf_space.n_dim) ').']);
    end
    
    %% metric checks    
    % Ensure presence of a metric that is a square symmetric matrix
    if ~isfield(s,'metric')
		s.metric = @(a,varargin) repmat(eye(size(shape_test_list,1)),size(a));
    else
        % test evaluation of the matrix
        metric_test = s.metric(shape_test_list{:});
        
        % size of test-evaluated matrix
        metric_size = size(metric_test);
        
        % Hard check for squareness
        if or(numel(metric_size) ~= 2, metric_size(1) ~= metric_size(2))
            error('Metric function must return a square matrix')
        % Hard check for appropriate size
        elseif metric_size(1) ~= n_shape
            error('Metric matrix is not of correct size for system')
        % Soft check for symmetry
        elseif ~(norm(metric_test-metric_test.',inf)<(norm(metric_test,inf)/10^6))
            warning('Metric matrix does not appear to be symmetric')
        end
    end
end