function output = sysf_three_link_lowRe(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Viscous Swimmer: 3-link'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Utilities/curvature_mode_toolbox/backbone_from_curvature_bases.m';...
				'Utilities/curvature_mode_toolbox/curvatures/discrete_joint_1.m';...
				'Utilities/curvature_mode_toolbox/curvatures/discrete_joint_2.m';...
				'Utilities/LowRE_toolbox/LowRE_dissipation_metric_from_curvature_bases.m'});

		case 'initialize'

			%Three-link Swimmer with 2:1 viscous drag ratio
            drag_ratio = 2;

            %Functional representations of local connection and metric

			s.A_num = @(alpha1,alpha2) Fullconnection3link(alpha1,alpha2,drag_ratio);
			s.metric =@(x,y) LowRE_dissipation_metric_from_curvature_bases...
				({@discrete_joint_1;@discrete_joint_2},[x;y],1,1,drag_ratio);

			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*2.5;

			%densities for various operations
			s.density.vector = [11 11 ]; %density to display vector field
			s.density.scalar = [21 21 ]; %density to display scalar functions
			s.density.eval = [21 21 ];   %density for function evaluations
            s.density.metric_eval = [11 11]; %density for metric evaluation
			s.finite_element_density = 21;

            %%%
			%Display parameters

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*1.5;
			s.tic_locs.y = [-1 0 1]*1.5;


			%%%%
			%Save the system properties
			output = s;


	end

end



function A=Fullconnection3link(alpha1,alpha2,drag_ratio)

% Specify geometry and fluid properties
l1=1/6; %Half-length of link 1
l2=1/6; %Half-length of link 2
l3=1/6; %Half-length of link 3

k=1;       %Fluid drag coefficient (only affects absolute cost of motion)

% Vector from base of link to midpoint
h1=[l1;0;0];
h2=[l2;0;0];
h3=[l3;0;0];

g_1 = groupmult(groupmult(-h2,[0,0,-alpha1]),-h1); %step back along middle link, rotate, then back along first link
gcirc_1_from_body = adjointinv2(g_1); %Jacobian from system body velocity to link 1 body velocity
gcirc_1_from_joint_1 = -adjointinv2(-h1); %Jacobian from joint 1 velocity to link 1 body velocity (negative sign is because of how alpha 1 is defined with respect to link 1
gcirc_1_from_joint_2 = zeros(3,3); %Jacobian from joint 1 velocity to link 1 body velocity
c1=[gcirc_1_from_body,gcirc_1_from_joint_1(:,3),gcirc_1_from_joint_2(:,3)]; %Total Jacobian from system velocities to link 1 body velocity

g_2 = eye(3);    % Link 2 is at the identity
gcirc_2_from_body = eye(3); %Jacobian from system body velocity to link 2 body velocity
gcirc_2_from_joint_1 = zeros(3,3); %Jacobian from joint 1 velocity to link 2 body velocity
gcirc_2_from_joint_2 = zeros(3,3); %Jacobian from joint 1 velocity to link 2 body velocity
c2=[gcirc_2_from_body,gcirc_2_from_joint_1(:,3),gcirc_2_from_joint_2(:,3)]; %Total Jacobian from system velocities to link 2 body velocity


g_3 = groupmult(groupmult(h2,[0,0,alpha2]),h3); %step forward along middle link, rotate, then forward along third link
gcirc_3_from_body = adjointinv2(g_3); %Jacobian from system body velocity to link 3 body velocity
gcirc_3_from_joint_1 = zeros(3,3); %Jacobian from joint 1 velocity to link 3 body velocity
gcirc_3_from_joint_2 = adjointinv2(h3); %Jacobian from joint 1 velocity to link 3 body velocity
c3=[gcirc_3_from_body,gcirc_3_from_joint_1(:,3),gcirc_3_from_joint_2(:,3)]; %Total Jacobian from system velocities to link 3 body velocity





forcemult1=[k*l1 0 0;
            0 drag_ratio*k*l1 0;
            0 0 k*(l1)^3];
        
forcemult2=[k*l2 0 0;
            0 drag_ratio*k*l2 0;
            0 0 k*(l2)^3];
        
forcemult3=[k*l3 0 0;
            0 drag_ratio*k*l3 0;
            0 0 k*(l3)^3];                                                                                                       
        
F=  [cos(alpha1), sin(alpha1), 0;
        -sin(alpha1), cos(alpha1), 0;
        l1*sin(alpha1), -l1*(1+cos(alpha1)), 1]*forcemult1*c1...
    ...
    +eye(3)*forcemult2*c2...
    ...
    +[cos(alpha2), -sin(alpha2), 0;
        sin(alpha2), cos(alpha2), 0;
        l3*sin(alpha2), l3*(1+cos(alpha2)), 1]*forcemult3*c3;

omega=F;
omega1=[omega(:,1),omega(:,2),omega(:,3)];
omega2=[omega(:,4),omega(:,5)];
A=omega1\omega2;
end

function a=adjointinv2(h)
a=[cos(h(3)) sin(h(3)) h(1)*sin(h(3))-h(2)*cos(h(3))
    -sin(h(3)) cos(h(3)) h(2)*sin(h(3))+h(1)*cos(h(3))
    0 0 1];
end

%SE(2) group multiplication of h*g for [x,y,theta] group elements
function a=groupmult(h,g)
a=[h(1)+g(1)*cos(h(3))-g(2)*sin(h(3));
    h(2)+g(2)*cos(h(3))+g(1)*sin(h(3));
    h(3)+g(3)];
end