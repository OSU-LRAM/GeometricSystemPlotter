function illustrate_backbone_shapespace(curvdef,paramvalues,orientation,plotnum)

% Specify orientation as midpoint-tangent unless specified otherwise
if ~exist('orientation','var')
    orientation = 'midpoint-tangent';
elseif or(isequal(orientation,'com-mean'), isequal(orientation,'midpoint-tangent'))
    % do nothing; orientation is a valid option
else
    % Get the user folder path
    load('sysplotter_config.mat','inputpath');
    % The hypothetical file path, supposing orientation is a system name
    calcfilePath = strcat(inputpath,'\sysplotter_data\sysf_',orientation,'_calc.mat');
    if exist(calcfilePath, 'file')
        % orientation = 'from-sys'; Orientation will be set to this within the
        % fatbackbone function.
    else
        warning(strcat('system data file sysf_',orientation,'_calc.mat can not be found. Defaulting to midpoint-tangent orientation...'))
        orientation = 'midpoint-tangent';
    end
end

% Identify what kind of backbone is being drawn
if isa(curvdef,'function_handle')
    fatbackbone_function = @fatbackbone_from_general_curvature;
elseif iscell(curvdef) && all(cellfun(@(x) isa(x,'function_handle'),curvdef))
    fatbackbone_function = @fatbackbone_from_curvature_bases;
else
    error('Unsupported curvature definition')
end

% Specify plot number as 171 unless specified otherwise
if ~exist('plotnum','var')
    plotnum = 171; % 171 was the hardcoded default before I modified this
end

% Get the gradient of parameter values in the first two dimensions
[~,grid_spacing_x] = gradient(paramvalues{1});
[grid_spacing_y,~] = gradient(paramvalues{2});

% Create the figure
fh = figure(plotnum);
close(fh);
fh = figure(plotnum);
set(gcf,'name',orientation);

axh = axes('Parent',fh);
axis equal
hold on 
% Set tick spacing
xticks(paramvalues{1}(:,1));
yticks(paramvalues{2}(1,:));
% Figure stylistic settings
box on % Repeat tickmarks and axis lines on the top and left sides

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