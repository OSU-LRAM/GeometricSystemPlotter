function output = shchf_kinsnake_demo_circle(input_mode,pathnames)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = 'Kinematic snake Circle Stroke';
			
		case 'dependency'
			
			output.dependency = {};
			
		case 'initialize'

			%%%%
			% Filename to save to

			%%
			%Path definitions

			%path definition
			p.phi_def{1} = {@strokedef};
			
			%marker locations
			p.phi_marker = [];
			
			%arrows to plot
			p.phi_arrows{1} = {2};

			%time to run path
			p.time_def{1} = {[0 2*pi]};

			p.cBVI_method{1}{1} = 'simple';

			%path resolution
			p.phi_res{1} = {1000};


			%%%%
			%Output the path properties
			output = p;
	end
	
end

function [stroke] = strokedef(t)

	t = t(:)';

	Rot=sqrt(2)/2*[1 1;-1 1];
	a=pi/8;

	stroke=(Rot*[-a*cos(t);-a*sin(t)])'+pi/4*[ones(size(t))', -ones(size(t))'];


end