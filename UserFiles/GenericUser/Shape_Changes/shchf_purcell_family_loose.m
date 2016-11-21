function output = shchf_purcell_family_loose(input_mode,pathnames)
 %[r1,r2] = convert.old_to_new_points(alpha1,alpha2);
	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
	
	switch input_mode
		
		case 'name'
			
			output = 'Family of square ''Purcell strokes''';
			
		case 'dependency'
			
			output.dependency = {};
			
		case 'initialize'

			%%%%
			%%
			%Path definitions
			
			%path definition
			A=linspace(0.5, 2, 5); % Amplitdude of circle
			
            
            % Multiple gaits (or piecewise definitions of gaits) can be
            % specified in a single file by making a nested
            % cell array for phi_def. The outer cell array contains
            % individual gaits, and the inner array contains segments of
            % the gaits. Here, we have multiple gaits with four segments
            % per gait
            %
            % Gait properties can be
            % specified for individual gaits by making the relevant fields
            % likewise cell structures, but a single value will be dealt
            % out to all the gaits and segments.
            
			for i = 1:numel(A)
				p.phi_def{i} = {...                 ith gait in family
                    @(t) strokedef1(t,A(i));...     first segment
                    @(t) strokedef2(t,A(i));...     second segment
                    @(t) strokedef3(t,A(i));...     third segment
                    @(t) strokedef4(t,A(i))};%      fourth segment		
				
				% Calculate the cBVI for this gait
%				p.cBVI_method{i}{1} = 'simple';
			end
			
			
			%marker locations
			p.phi_marker = [];
			
			%arrows to plot
			p.phi_arrows = 0;

			%time to run path
			p.time_def = [0 1];


			%path resolution
			p.phi_res = 50;


			%%%%
			%Output the shch properties
			output = p;

	end
	
end


% Segments are each one side of a square.
function [stroke] = strokedef1(t,A)

	t = t(:);
	
	stroke = A*[ones(size(t)), (1-2*t)];
	
	
end

function [stroke] = strokedef2(t,A)

	t = t(:);
	
	stroke = A*[(1-2*t), -ones(size(t))];
	
	
end

function [stroke] = strokedef3(t,A)

	t = t(:);
	
	stroke = A*[-ones(size(t)), -(1-2*t)];
	
	
end

function [stroke] = strokedef4(t,A)

	t = t(:);
	
	stroke = A*[-(1-2*t), ones(size(t))];
	
	
end