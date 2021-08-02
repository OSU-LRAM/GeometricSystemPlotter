function output = shchf_circle_family_gen(input_mode,pathnames)
	% Default argument
	if ~exist('input_mode','var')
		input_mode = 'initialize';
    end
	switch input_mode
		case 'name'
			output = 'General-case family of circular strokes for three-link swimmers';
		case 'dependency'
			output.dependency = {};
		case 'initialize'
			%path definition
			A=0.25:0.25:2; %amplitude of circle
            phi=pi/4:pi/4:2*pi; %starting phase
			
            % define gait family
            ctr = 0;
			for i = 1:length(A)
                for j = 1:length(phi)
                    ctr = ctr + 1;
                    p.phi_def{ctr,1}{1} = @(t) strokedef(t,A(i),phi(j));		

                    % Calculate the cBVI for this gait
                    p.cBVI_method{ctr}{1} = 'simple';
                end
            end
			
			%marker locations
			p.phi_marker = [];
			
			%arrows to plot
			p.phi_arrows = 1; %[repmat({{0}},numel(A)-1,1);{{0}}];

			%time to run path
			p.time_def = [0 2*pi];

			%path resolution
			p.phi_res = 70;

			%%%%
			%Output the shch properties
			output = p;
    end
end

function [stroke] = strokedef(t,A,phi)
	t = -t(:)'; %follow gait CCW
	stroke = (sqrt(2)/2*[1 -1;1 1]*(A*[cos(t+phi); sin(t+phi)]))';
end