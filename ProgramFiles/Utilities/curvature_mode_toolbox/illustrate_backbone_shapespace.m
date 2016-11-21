function illustrate_backbone_shapespace(curvdef,paramvalues)

% Get the gradient of parameter values in the first two dimensions
[~,grid_spacing_x] = gradient(paramvalues{1});
[grid_spacing_y,~] = gradient(paramvalues{2});

fh = figure(171);
close(fh);
fh = figure(171);


axh = axes('Parent',fh);
axis equal
hold on

for idx = 1:numel(paramvalues{1})
    
    % Extract the parameter values for this combination of parameters
    p = zeros(size(paramvalues));
    for idx2 = 1:numel(paramvalues)
        p(idx2) = paramvalues{idx2}(idx);
    end
    
    % The backbone should be as long as the spacing between the nearby
    % elements
    blnth = .9 * max(grid_spacing_x(idx),grid_spacing_y(idx));
    
    %Generate the backbone locus
    B = fatbackbone_from_general_curvature(curvdef,p,blnth,blnth/20);
    
    for idx2 = 1:2
        B(:,idx2) = B(:,idx2) + p(idx2);% - mean(B(:,idx2));
    end
    
    % draw the backbone at the specified location
    plot(B(:,1),B(:,2),'Parent',axh)

end