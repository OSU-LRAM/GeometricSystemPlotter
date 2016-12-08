%If metric is defined with separate numerator and denominator, merge
%them. 
function s = merge_metric(s)

    if isfield(s,'metric_den')

        %get the list of connection vector field zoom levels
        zoom_list = fieldnames(s.metricfield);

        %Get the size of the local connection matrix
        M_size = size(s.metricfield.(zoom_list{1}).content.metric);

        %loop over all zoom levels of the field being evacuated
        for m = 1:length(zoom_list)

            % Allocate a cell array to hold the local connection
            s.metricfield.(zoom_list{m}).content.metric = cell(M_size);

            %first, apply the smart divider to
            %get the base connection vector field Avec and the location of any
            %singularities
            [s.metricfield.(zoom_list{m}).content.metric ...
                ,s.metricfield.(zoom_list{m}).singularities]...
                = cellfun(@(num,den) smart_divider(s.grid.(s.metricfield.(zoom_list{m}).type),num,den)...
                ,s.metricfield.(zoom_list{m}).content.metric...
                ,s.metricfield.(zoom_list{m}).content.metric_den, 'UniformOutput',false);

        end
        
    end

end