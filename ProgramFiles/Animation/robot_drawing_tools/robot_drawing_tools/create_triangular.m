function robot = create_triangular(parent,data_source,system_name,orientation_color)

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



%Get a serpenoid body
% body_local = fatbackbone_from_curvature_bases({@serpenoid_1;@serpenoid_2},[0;0],1,.2/6);
body_local = fatbackbone_from_curvature_bases_triangular({@serpenoid_1;@serpenoid_2},[0;0],1,.2/6);
% body_local1 = backbone_from_curvature_bases({@serpenoid_1;@serpenoid_2},[0;0],1);

%%%
%build the robot

%%
%Make the polygons for the links (really just one curved link for drawing
%purposes)
robot.links = object_primative;

%Size of links
link_scale = 1;

for i = 1:size(robot.links,1)
	robot.links(i).type = 'polygon';
	robot.links(i).local_geom.vertices =  link_scale * body_local;
	robot.links(i).graphics.fill = {'FaceColor','w'};
	robot.links(i).graphics.edge = {'EdgeColor','k','LineWidth',1};
end



%make the center-dot
robot.center = line('XData',[],'YData',[],'Color','k','LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor','k');

%make the orientation line
if ~exist('orientation_color','var')
	orientation_color = [110 110 110]/255;
end
robot.orientation = line('XData',[],'YData',[],'Color',orientation_color,'LineStyle','-','Marker','none','LineWidth',2);

% Identify the system file containing optimal coordinate information
robot.sysdata = fullfile(data_source,['sysf_' system_name '_calc']);