function illustrate_shapespace(system,plotnum)

geometry = system.geometry;
visual = system.visual;

load('sysplotter_config','Colorset')

        
% Specify plot number as 177 unless specified otherwise
if ~exist('plotnum','var')
    plotnum = 177;
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
set(fh,'name',['Baseframe: ' geometry.baseframe]);

axh = axes('Parent',fh);
axis(axh,'equal')
hold(axh,'on')
% Set tick spacing
axh.XTick = paramgrid{1}(:,1);
axh.YTick = paramgrid{2}(1,:);
set(axh,'XLim', [min(paramgrid{1}(:))-grid_spacing_x(1)/2,max(paramgrid{1}(:))+grid_spacing_x(end)/2]);
set(axh,'YLim', [min(paramgrid{2}(:))-grid_spacing_y(1)/2,max(paramgrid{2}(:))+grid_spacing_y(end)/2]);
% Figure stylistic settings
axh.Box = 'on'; % Repeat tickmarks and axis lines on the top and left sides
grid(axh,'on'); % show a grid at the specified shape points

% The system length scale should be as long as the minimum spacing between nearby
% elements
blnth = .9 * min([grid_spacing_x(:);grid_spacing_y(:)]);

% Scale the geometry defined in the system file
 geometry.length = blnth;


% Specfication of face colors for different parts of the system. Currently
% support primary portion of the system in white with a black outline, and
% secondary portion in the spot color with a black outline, and a third
% portion in the secondary color
fillcolors = {'w',Colorset.spot,Colorset.secondary};

 
for idx = 1:numel(paramgrid{1})
    
    % Extract the parameter values for this combination of parameters
    p = zeros(size(paramgrid));
    for idx2 = 1:numel(paramgrid)
        p(idx2) = paramgrid{idx2}(idx);
    end
    
    
    %Generate the backbone locus
    B = generate_locomotor_locus(geometry,p,visual);
    
    % Force B to be a column cell array
    if ~iscell(B)
        B = {B};
    end
    B = B(:);
    
    for idx2 = 1:numel(B)
        
        % Iterate over each piece within the body
        
        % Force B{idx2} to be a column cell array
        if ~iscell(B{idx2})
            B{idx2} = {B{idx2}};
        end
        B{idx2} = B{idx2}(:);
        
        for idx3 = 1:numel(B{idx2})
            % transform the body elements to their position on the grid
            for idx4 = 1:2
                B{idx2}{idx3}(idx4,:) = B{idx2}{idx3}(idx4,:) + p(idx4);
            end
    
            % draw the backbone at the specified location
            bp{idx,idx2,idx3} = patch('Xdata',B{idx2}{idx3}(1,:),'Ydata',B{idx2}{idx3}(2,:),'Parent',axh,'EdgeColor','k','FaceColor',fillcolors{idx2});
        end
    end
end