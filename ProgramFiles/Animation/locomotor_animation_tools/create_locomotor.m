function robot = create_locomotor(...
    parent,...
    data_source,...
    system_name,...
    info_needed)


% Extract geometry and visualization information from system
s = info_needed.s;
geometry = s.geometry;
n_shape = nargin(s.A_num);

if isfield(s,'visual')
    visual = s.visual;
else
    visual = struct();
end

% Get the system colors
load('sysplotter_config','Colorset');

%Basic object
object_primative = struct('parent',parent,...
	'type',[],...
    'local_geom',struct('vertices',[],'scaling',[],'sidelengths',[],'radius',[],'corner_rad',[]),...
    'ref_pos', struct('now',[],...
        'last',[]...
        ),...
    'geom',struct('vertices',[]...
        ,'line_eq',struct('matrix',[])...
        ,'segments',[]...
        ),...
    'graphics',struct('handle',[],'fill',[],'edge',[])...
    );

%%%%%%%
% Generate the robot geometry at a test shape value

% Choose a random point within the plot range (to probabalistically
% avoid testing at a singularity)
shape_test_point = ones(n_shape,1);
for idx = 1:n_shape
   shape_test_point(idx) = ...
       (rand(1)*(s.grid_range(2*idx)-s.grid_range((2*idx)-1)))-s.grid_range((2*idx)-1);
end

% Get a locus of points for the shapes with which to draw the locomotor.
% This may be either an array of points, or a cell array of point arrays,
% each of which would be drawn separately
B = generate_locomotor_locus(geometry,shape_test_point,visual);

% Force B to be a column cell array with a cell inside it
if ~iscell(B)
    B = {B};
end
B = B(:);

if ~iscell(B{1})
    B{1} = {B{1}};
end
B{1} = B{1}(:);

%%%
% Make a cell array in which each cell is one element of the drawing that
% should be made (taken from the dimensionality of B. A simple backbone or fat chain would be a single cell
% and multiple cells could be used, for example, to draw wheels next to the
% link in a different visual style from the body
robot.body = repmat(object_primative,size(B));

%%%
% Specify how to draw the locomotor

% Specfication of face colors for different parts of the system. Currently
% support primary portion of the system in white with a black outline, and
% secondary portion in the spot color with a black outline, and a third
% portion in the secondary color
fillcolors = {'w',Colorset.spot,Colorset.secondary};

for i = 1:size(robot.body,1)
	robot.body(i).type = 'polygon';
	robot.body(i).local_geom.vertices =  B{i};
	robot.body(i).graphics.fill = {'FaceColor',fillcolors{i}};
	robot.body(i).graphics.edge = {'EdgeColor','k','LineWidth',1};
end



%make the center-dot
robot.center = line(...
    'XData',[],...
    'YData',[],...
    'Color','k',...
    'LineStyle','none',...
    'Marker','o',...
    'MarkerSize',10,...
    'MarkerFaceColor',Colorset.spot,...
    'parent', parent);

robot.orientation = line(...
    'XData',[],...
    'YData',[],...
    'Color',Colorset.secondary,...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'parent', parent);

% Identify the system file containing optimal coordinate information
robot.sysdata = fullfile(data_source,['sysf_' system_name '_calc']);