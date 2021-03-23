function output = shchf_3d_ellipse(input_mode,pathnames)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = '3D Ellipse, Tunable';
			
		case 'dependency'
			
			output.dependency = {};
			
		case 'initialize'

			%%%%
			% Filename to save to

			%%
			%Path definitions

			%path definition
			p.phi_def = @strokedef;
			
			%marker locations
			p.phi_marker = [];
			
			%arrows to plot
			p.phi_arrows = 2;

			%time to run path
			p.time_def = [0 2*pi];

			%p.cBVI_method{1}{1} = 'simple';

			%path resolution
			p.phi_res = 100;


			%%%%
			%Output the path properties
			output = p;
	end
	
end

function [stroke] = strokedef(t)

	t = -t(:);

    yaw = 80 * (pi/180);
    pitch = 20 * (pi/180);
    roll = 55 * (pi/180);
    
    x_o = 0;
    y_o = 0;
    z_o = .2;
    
    yaw_R = [cos(yaw),sin(yaw),0,0;-sin(yaw),cos(yaw),0,0;0,0,1,0;0,0,0,1];
    pitch_R = [1,0,0,0;0,cos(pitch),sin(pitch),0;0,-sin(pitch),cos(pitch),0;0,0,0,1];
    roll_R = [cos(roll),0,-sin(roll),0;0,1,0,0;sin(roll),0,cos(roll),0;0,0,0,1];
    
	a=12;
    b=7;

	stroke=[-a*cos(t)+x_o,-b*sin(t)+y_o,zeros(numel(t),1)+z_o,ones(numel(t),1)]*yaw_R*pitch_R*roll_R;
    stroke = stroke(:,1:3);


end