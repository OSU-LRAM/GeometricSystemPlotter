function illustrate_backbone_shapespace(curvdef,paramvalues,orientation)

% Specify orientation as midpoint-tangent unless specified otherwise
if ~exist('orientation','var')
    orientation = 'midpoint-tangent';
end

% Identify what kind of backbone is being drawn
if isa(curvdef,'function_handle')
    fatbackbone_function = @fatbackbone_from_general_curvature;
elseif iscell(curvdef) && all(cellfun(@(x) isa(x,'function_handle'),curvdef))
    fatbackbone_function = @fatbackbone_from_curvature_bases;
else
    error('Unsupported curvature definition')
end

%Make sure that paramvalues are in a column cell to avoid any row/column
%problems below
paramvalues = paramvalues(:);

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
    B = fatbackbone_function(curvdef,p,blnth,blnth/20,orientation);
    
    for idx2 = 1:2
        B(:,idx2) = B(:,idx2) + p(idx2);
    end
    
    % draw the backbone at the specified location
    plot(B(:,1),B(:,2),'Parent',axh,'Color','k')

end