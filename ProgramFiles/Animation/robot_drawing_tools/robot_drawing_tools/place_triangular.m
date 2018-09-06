function robot = place_triangular(robot,config,coords)
%set the positions of the links and wheelsets of the robot in a given
%configuration and coordinate scheme

alpha_1 = config.alpha_1;
alpha_2 = config.alpha_2;
x = config.x;
y = config.y;
theta = config.theta;


% Rebuild the backbone
% body_local = fatbackbone_from_curvature_bases({@serpenoid_1;@serpenoid_2},[alpha_1;alpha_2],1,.2/6);
body_local = fatbackbone_from_curvature_bases_triangular({@serpenoid_1;@serpenoid_2},[alpha_1;alpha_2],1,.2/6);
robot.links.local_geom.vertices =  body_local;

		
%Build the transform for the position
ref_pos = vec_to_mat_se2([x y theta]);

robot.ref_pos = [x y theta];

%Build the transform from the body frame to the middle link frame
switch coords
	
	case 'original'
		
		coord_pos = eye(3);
		
	case 'mean'
		
		beta_theta = (alpha_2-alpha_1)/3;
		
		coord_pos = vec_to_mat_se2([0 0 -beta_theta]);
		
	case 'optimized'
		
		load(robot.sysdata,'s')
		
		beta_x = interpn(s.grid.eval{1}, s.grid.eval{2}, s.B_optimized.eval.Beta{1},alpha_1,alpha_2);
		beta_y = interpn(s.grid.eval{1}, s.grid.eval{2}, s.B_optimized.eval.Beta{2},alpha_1,alpha_2);
		beta_theta = interpn(s.grid.eval{1}, s.grid.eval{2}, s.B_optimized.eval.Beta{3},alpha_1,alpha_2);
		
		coord_pos = inv(vec_to_mat_se2([beta_x, beta_y, beta_theta]));
		
	case 'original-to-mean'
		
		beta_theta = (alpha_2-alpha_1)/3*config.percentage;
		
		coord_pos = vec_to_mat_se2([0 0 -beta_theta]);
		
	otherwise
		
		error('Unsupported coordinate set');
		
end

%%%%%%
%Get the world-positions of the robot elements
	
linkposmat = ref_pos*coord_pos;


%Pull the position vectors back out
linkpos = mat_to_vec_se2(linkposmat)';

%Assign the link positions to the link and wheelset polygons
for i = 1:size(linkpos,2)
	
	robot.links(i).ref_pos.now = deal(linkpos(:,i));

end

% %Put the center-dot at the ref position
% set(robot.center,'XData',x,'YData',y)
