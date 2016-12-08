function output = shchf_LowRe_MaxEff_and_Disp(input_mode,pathnames)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	switch input_mode
		
		case 'name'
			
			output = 'Low Re max eff and disp strokes';
			
		case 'dependency'
			
			output.dependency = {}; % No dependencies
			
		case 'initialize'

			%%%%
			% Filename to save to

			%%
			%Path definitions

            % Multiple gaits (or piecewise definitions of gaits) can be
            % specified in a single file by making a nested
            % cell array for phi_def. The outer cell array contains
            % individual gaits, and the inner array contains segments of
            % the gaits. Here, we have two gaits but only one segment
            % per gait.
            %
            % Gait properties can be
            % specified for individual gaits by making the relevant fields
            % likewise cell structures, but a single value will be dealt
            % out to all the gaits and segments.
			p.phi_def{1} = {@strokedef_eff};
			p.phi_def{2} = {@(t)strokedef_disp(t+pi/2)};
			
			%marker locations
			p.phi_marker = [];
			
			%arrows to plot
			p.phi_arrows = {{0} {2}};

			%time to run path
			p.time_def = [0 2*pi];


			%path resolution
			p.phi_res = 50;


			%%%%
			%Save the shape change properties
			output = p;

	end
	
end

function [stroke] = strokedef_eff(t)

	t = -t(:)';

	XV=[ 1.495398722861924 1.326600886340017 0.083253854456986 0.089227148110158  -0.034695750579259  -0.044759551791058];
	Rot=sqrt(2)/2*[1 -1;1 1];
	a1=XV(1); b1=XV(2); a3=XV(3); b3=XV(4); a5=XV(5); b5=XV(6);

	stroke=(Rot*[a1*cos(t)+a3*cos(3*t)+a5*cos(5*t);b1*sin(t)+b3*sin(3*t)+b5*sin(5*t)])';


end

function [stroke] = strokedef_disp(t)

	t = -t(:)';

XV=[2.017955906441066   1.353385144748382   0.013112519773635   0.263823759269392  -0.031456173657342  -0.003626287235861];	Rot=sqrt(2)/2*[1 -1;1 1];
	a1=XV(1); b1=XV(2); a3=XV(3); b3=XV(4); a5=XV(5); b5=XV(6);

	stroke=(Rot*[a1*cos(t)+a3*cos(3*t)+a5*cos(5*t);b1*sin(t)+b3*sin(3*t)+b5*sin(5*t)])';


end