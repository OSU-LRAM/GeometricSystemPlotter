function robot = place_locomotor(robot,config,s)
% Use the locomotor geometry and its current configuration to place its
% body locus in the world

% Extract the configuration parameters
shapeparams = config.shape;
x = config.x;
y = config.y;
theta = config.theta;

% Get the geometry and visualization parameters for the system
geometry = s.geometry;
if isfield(s,'visual')
    visual = s.visual;
else
    visual = struct();
end


%%%%%%
% Rebuild the backbone

% Evaluate points on the locomotor locus
B = generate_locomotor_locus(geometry,shapeparams,visual);

% Force B to be a column cell array
if ~iscell(B)
    B = {B};
end
B = B(:);

% Load the points into the robot structure
for i = 1:size(robot.body,1)
	robot.body(i).local_geom.vertices =  B{i}';
end

%Build the transform for the position
robot.ref_pos = [x y theta];

%Assign this position to each class of items in the body
for i = 1:size(robot.body)
	
	robot.body(i).ref_pos.now = robot.ref_pos;%deal(linkpos(:,i));

end