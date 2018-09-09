function illustrate_shapespace(system,plotnum)

geometry = system.geometry;
visual = system.visual;


% Identify what kind of system is being drawn, and use this to specify how
% the drawing should be generated
switch geometry.type
    
    case {'curvature basis','curvature bases','general curvature'}
        drawing_generator = @fat_backbone;
        
    case {'n-link chain'}
        drawing_generator = @fat_chain;
        
end
        
% Specify plot number as 171 unless specified otherwise
if ~exist('plotnum','var')
    plotnum = 171;
end

% Take the display grid from the provided display structure. Make sure that
% the components of the grid are in a column cell array to avoid any
% row/column problems below
paramgrid = visual.grid(:);

% Get the spacing of parameter values in the first two dimensions
[~,grid_spacing_x] = gradient(paramgrid{1});
[grid_spacing_y,~] = gradient(paramgrid{2});

% Create the figure
fh = figure(plotnum);
clf(fh);
fh = figure(plotnum);
set(fh,'name',geometry.baseframe);

axh = axes('Parent',fh);
axis(axh,'equal')
hold(axh,'on')
% Set tick spacing
axh.XTick = paramgrid{1}(:,1);
axh.YTick = paramgrid{2}(1,:);
% Figure stylistic settings
axh.Box = 'on'; % Repeat tickmarks and axis lines on the top and left sides
grid(axh,'on'); % show a grid at the specified shape points

% The system length scale should be as long as the minimum spacing between nearby
% elements
blnth = .9 * min([grid_spacing_x(:);grid_spacing_y(:)]);

% Scale the geometry defined in the system file
 geometry.length = blnth;


for idx = 1:numel(paramgrid{1})
    
    % Extract the parameter values for this combination of parameters
    p = zeros(size(paramgrid));
    for idx2 = 1:numel(paramgrid)
        p(idx2) = paramgrid{idx2}(idx);
    end
    
    
    %Generate the backbone locus
    B = drawing_generator(geometry,p,visual);
    
    for idx2 = 1:2
        B(idx2,:) = B(idx2,:) + p(idx2);
    end
    
    % draw the backbone at the specified location
    plot(B(1,:),B(2,:),'Parent',axh,'Color','k')
    

    
end