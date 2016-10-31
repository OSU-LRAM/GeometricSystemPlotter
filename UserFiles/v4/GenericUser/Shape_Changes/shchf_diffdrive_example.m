function output = shchf_diffdrive_example(input_mode)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = 'Diffdrive Example Motion';
			
		case 'dependency'
			
			output.dependency = {};
			
		case 'initialize'

			%%
			%Path definitions

			%path definition. Concatenate these two shape change segments into a
			%single trajectory (best results when the end of one segment
			%starts the next.
			p.phi_def{1} = {@drive_forward, @turn_in_place}; 
			
			
			%marker locations
			p.phi_marker = []; % No marker on this path (can put, e.g. endpoints of path if desired)
			
			%arrows to plot
			p.phi_arrows{1} = {1,2}; %put one direction arrow on the first segment, two on the second

			%time to run path
			p.time_def{1} = {[0 1], [0 1]}; % Duration of each segment

			% With a closed-loop gait, enabling this line takes the area
			% integral of the Constraint Curvature Function (note that this
			% is currently a slow operation.
			%p.cBVI_method{1}{1} = 'simple';

			%number of points in each path.
			p.phi_res{1} = {100, 100};


			%%%%
			%Output the path properties
			output = p;
	end
	
end

function [alpha] = drive_forward(t)

	t = t(:);

	alpha = [-1+t -1+t];


end


function [alpha] = turn_in_place(t)

	t = t(:);

	alpha = [-t t];


end