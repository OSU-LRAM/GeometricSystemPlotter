function [metric, metricgrad] = metric_interpolator(y,s,n,dimension)
%% Preliminaries for gradient calculation
% Preallocating memory for variables which we will need in further
% calculation 
yvalues=cell(n,dimension); % Cell representation of the coordinates of all points forming the gait
interpmetricgrid=cell(1,dimension);  % Variable which will store the metric grid used for interpolation
metric1=zeros(n,dimension,dimension); % Variable which will store metric at each point in the form of a matrix
metric = repmat({zeros(dimension)},[n 1]); % Variable which stores the metric at each point in the form of a 2x2 matrix
metricgrad1=zeros(n,dimension,dimension,dimension); % Variable which will store gradient of metric at each point in the form of a matrix
metricgrad = repmat({zeros(dimension)},[n,dimension]); % Variable which will store gradient of metric at each point in the form of a matrix
afactor=0.001;

% Interpolation to calculate all the variables needed for gradient
% calculation
for i=1:1:n
    for j=1:1:dimension
        yvalues{i,j}=y(i,j);
    end
end

y_for_interp = mat2cell(y,size(y,1),ones(1,size(y,2)));

for j=1:1:dimension
    interpmetricgrid{j}=s.grid.metric_eval{j,1};
end

for j=1:1:dimension
    for k=1:1:dimension
        metric1(:,j,k)=interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y_for_interp{:},'spline');
    end
end

if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'acceleration coord')
    for i=1:n
       metric{i}=eye(dimension);
    end
elseif strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
    for i=1:n
        for j=1:1:dimension
           for k=1:1:dimension
               metric{i}(j,k)=metric1(i,j,k);
           end
        end
    end
end
if strcmpi(s.costfunction,'pathlength metric2')
    for i = 1:n
        metric{i} = metric{i}*metric{i};
    end
end

if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'acceleration coord')
    for i=1:dimension
        for j=1:n
            metricgrad{j,i} = zeros(dimension);
        end
    end
elseif strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
    y2 = zeros(size(y));
    y1 = y2;
    for l=1:1:dimension
        for m=1:1:dimension
            if m==l
               y2(:,m)=y(:,m)+afactor*ones(length(y),1);
               y1(:,m)=y(:,m)-afactor*ones(length(y),1);
            else
               y2(:,m)=y(:,m);
               y1(:,m)=y(:,m);
            end
        end
        y2_for_interp = mat2cell(y2,size(y,1),ones(1,size(y,2)));
        y1_for_interp = mat2cell(y1,size(y,1),ones(1,size(y,2)));
        for j=1:1:dimension
            for k=1:1:dimension
                metricgrad1(:,l,j,k)=(interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y2_for_interp{:},'spline')...
                    -interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y1_for_interp{:},'spline'))/(2*afactor);
            end
        end
        for i=1:n
            for j=1:1:dimension
                for k=1:1:dimension
                    metricgrad{i,l}(j,k)=metricgrad1(i,l,j,k);
                end
            end
        end
    end
end
if strcmpi(s.costfunction,'pathlength metric2')
    for l = 1:dimension
        for i = 1:n
            metricgrad{i,l} = metricgrad{i,l}*metric{i}+metric{i}*metricgrad{i,l};
        end
    end
end


end