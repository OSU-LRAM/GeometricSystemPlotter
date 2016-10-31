%Evaluate the metric over the fine grid for calculations and the coarse
%grid for vector display
function s = evaluate_metric(s)

	% Get the number of shape variables
	n_shape = nargin(s.A_num);
	shape_test_list = num2cell(ones(n_shape,1));


	%Ensure presence of a metric and a metric denominator
	if ~isfield(s,'metric')
		s.metric = @(varargin) eye(size(shape_test_list,1));
	end

	if ~isfield(s,'metric_den')
		s.metric_den = @(varargin) ones(size(shape_test_list,1));
	end

    %list of zoom levels at which to evaluate connection vector fields and the
    %grid which should be used for them
    metric_field_list = {'display','vector';
                         'metric_eval','metric_eval'};
                     
    %loop over list, creating the vector fields
    for i = 1:size(metric_field_list,1);
        
        %generate one large array with all the metric information in it for
        %each of the numerator and denominator, and treat the set as two
        %cells in an array
        %note that the fields are the _negative_ of the connection
        s.metricfield.(metric_field_list{i,1}).content.metric_num =... %numerator
            arrayfun(@(varargin) s.metric(varargin{:}),s.grid.(metric_field_list{i,2}){:},'UniformOutput',false);
        s.metricfield.(metric_field_list{i,1}).content.metric_den =... %denominator
            arrayfun(@(varargin) s.metric_den(varargin{:}),s.grid.(metric_field_list{i,2}){:},'UniformOutput',false);
        
        %mark what grid was used to create the field
        s.metricfield.(metric_field_list{i,1}).type = (metric_field_list{i,2});
        
    end

end