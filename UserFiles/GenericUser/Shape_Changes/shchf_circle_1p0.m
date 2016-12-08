function output = shchf_circle_1p0(input_mode,pathnames)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = 'Circle Stroke, 1.0 amplitude';
			
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

			p.cBVI_method = 'simple';

			%path resolution
			p.phi_res = 100;


			%%%%
			%Output the path properties
			output = p;
	end
	
end

function [stroke] = strokedef(t)

	t = -t(:)';

	Rot=sqrt(2)/2*[1 -1;1 1];
	a=1.0;

	stroke=(Rot*[-a*cos(t);-a*sin(t)])';


end