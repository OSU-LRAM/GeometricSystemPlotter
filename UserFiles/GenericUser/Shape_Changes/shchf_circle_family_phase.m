function output = shchf_circle_family_phase(input_mode,pathnames)


	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
	
	switch input_mode
		
		case 'name'
			
			output = 'Family of circular strokes (by phase) for three-link swimmers';
			
		case 'dependency'
			
			output.dependency = {};
			
		case 'initialize'
			
			%path definition
			P= 0:pi/2:2*pi; % starting phase
			
            % Multiple gaits (or piecewise definitions of gaits) can be
            % specified in a single file by making a nested
            % cell array for phi_def. The outer cell array contains
            % individual gaits, and the inner array contains segments of
            % the gaits. Here, we have multiple gaits but only one segment
            % per gait.
            %
            % Gait properties can be
            % specified for individual gaits by making the relevant fields
            % likewise cell structures, but a single value will be dealt
            % out to all the gaits and segments.
			for i = 1:1:numel(P)
				p.phi_def{i,1}{1} = @(t) strokedef(t,P(i));		
				
				% Calculate the cBVI for this gait
				p.cBVI_method{i}{1} = 'simple';
			end
			
			
			%marker locations
			p.phi_marker = [];
			
			%arrows to plot
			p.phi_arrows = 0; %[repmat({{0}},numel(A)-1,1);{{0}}];

			%time to run path
			p.time_def = [0 2*pi];


			%path resolution
			p.phi_res = 70;


			%%%%
			%Output the shch properties
			output = p;

	end
	
end

function [stroke] = strokedef(t,P)

	t = t(:)';

	stroke = (sqrt(2)/2*[1 -1;1 1]*([cos(t + P); sin(t + P)]))';

end