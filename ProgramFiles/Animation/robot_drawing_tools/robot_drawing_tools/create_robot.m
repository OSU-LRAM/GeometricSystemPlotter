function robot = create_robot(parent,data_source,system_name,orientation_color)

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

%Size of links
link_scale = 1/6/1.05;

%relative scale of wheels to links
wheel_scale = .3;

%placement of wheels relative to links
wheel_outboard = .2*[1;-1]*link_scale;

%Get a 50-point squashed rectangle for body elements
[X_body,Y_body] = squashed_rectangle(2,.2,.3,50);

%Get a 50-point squashed rectangle for wheel elements
[X_wheel,Y_wheel] = squashed_rectangle(2*wheel_scale,.4*wheel_scale,.2,50);

%%%
%build the robot

%%
%Make the polygons for the links
robot.links = [object_primative;object_primative;object_primative];



for i = 1:size(robot.links,1)
	robot.links(i).type = 'polygon';
	robot.links(i).local_geom.vertices =  link_scale * [X_body Y_body];
	robot.links(i).graphics.fill = {'FaceColor','w'};
	robot.links(i).graphics.edge = {'EdgeColor','k','LineWidth',1};
end

%%
%Make the polygons for the wheelsets
robot.wheels = [object_primative object_primative;...
	object_primative object_primative;...
	object_primative object_primative];


for i = 1:size(robot.wheels,1)
	
	for j = 1:size(robot.wheels,2)
		
		robot.wheels(i,j).type = 'polygon';
		robot.wheels(i,j).local_geom.vertices =  link_scale*[X_wheel Y_wheel];
		robot.wheels(i,j).local_geom.vertices(:,2) = robot.wheels(i,j).local_geom.vertices(:,2) + wheel_outboard(j);
		robot.wheels(i,j).graphics.fill = {'FaceColor',[234 14 30]/255};
		robot.wheels(i,j).graphics.edge = {'EdgeColor','k'};
		
	end
	
end

%make the center-dot
robot.center = line('XData',[],'YData',[],'Color','k','LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor','k');

if exist('orientation_color')
	%make the orientation line
	robot.orientation = line('XData',[],'YData',[],'Color',orientation_color,'LineStyle','-','Marker','none','LineWidth',2);
end

% Identify the system file containing optimal coordinate information
robot.sysdata = fullfile(data_source,['sysf_' system_name '_calc']);