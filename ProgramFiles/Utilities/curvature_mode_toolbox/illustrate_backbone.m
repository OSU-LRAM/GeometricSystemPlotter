function illustrate_backbone(curvdef,paramvalue,orientation,plotnum)
% taken from illustrate_backbone_shapespace
% plots a single snake

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

% Create the figure
fh = figure(plotnum);
set(gcf,'name',orientation);

axh = axes('Parent',fh);
axis equal
axh.XLim = [-1 1];
axh.YLim = [-1 1];
hold on 

% Figure stylistic settings
box on % Repeat tickmarks and axis lines on the top and left sides

%for idx = 1:numel(paramvalues{1})
idx=1;
    
    p=paramvalue;
    
    % The backbone should be as long as the spacing between the nearby
    % elements
    blnth = 1;
    
    %Generate the backbone locus
    B = fatbackbone_function(curvdef,p,blnth,blnth/20,orientation);
    
    % draw the backbone at the specified location
    plot(B(:,1),B(:,2),'Parent',axh,'Color','k')
