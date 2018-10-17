function output = shchf_LowRe_MaxEfficiency(input_mode,pathnames)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = 'Low Re maximum-effiency stroke';
			
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

	XV=[ 1.495398722861924 1.326600886340017 0.083253854456986 0.089227148110158  -0.034695750579259  -0.044759551791058];
	Rot=sqrt(2)/2*[1 -1;1 1];
	a1=XV(1); b1=XV(2); a3=XV(3); b3=XV(4); a5=XV(5); b5=XV(6);

	stroke=(Rot*[a1*cos(t)+a3*cos(3*t)+a5*cos(5*t);b1*sin(t)+b3*sin(3*t)+b5*sin(5*t)])';


end