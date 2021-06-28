function output = shchf_diffdrive_cyclic_aligned(input_mode,pathnames)
	% Default argument
    if ~exist('input_mode','var')
		input_mode = 'initialize';
    end
    
    switch input_mode
		case 'name'
			output = 'Diffdrive Cyclic Motion, Aligned';
		case 'dependency'
			output.dependency = {};
		case 'initialize'
            % path definition
            % equivalent to example motion, but cyclic (parallel parking)
            side_lengths = pi/32:pi/64:pi/8;
            p.phi_def = cell(1,length(side_lengths));
            for i = 1:length(side_lengths)
                p.phi_def{i} = {@(t) drive_forward(t, side_lengths(i)),...
                                @(t) turn_left(t, side_lengths(i)),...
                                @(t) drive_backward(t, side_lengths(i)),...
                                @(t) turn_right(t, side_lengths(i))};
            end
            
			%marker locations
			p.phi_marker = []; % No marker on this path (can put, e.g. endpoints of path if desired)
			
			%arrows to plot
			p.phi_arrows = {{1,1,1,1}};

			%time to run path
			p.time_def = [0 1]; % Duration of each segment

			% enable area integration
			p.cBVI_method{1}{1} = 'simple';

			%number of points in each path.
			p.phi_res = 20;

			% output the path properties
			output = p;
    end
end

function [alpha] = drive_forward(t,len)
	t = t(:)*len;
    start = [0 0];
	alpha = start + [t zeros(size(t))];
end

function [alpha] = turn_left(t,len)
	t = t(:)*len;
	alpha = drive_forward(1,len) + [zeros(size(t)) t];
end

function [alpha] = drive_backward(t,len)
    t = t(:)*len;
    alpha = turn_left(1,len) + [-t zeros(size(t))];
end

function [alpha] = turn_right(t,len)
    t = t(:)*len;
    alpha = drive_backward(1,len) + [zeros(size(t)) -t];
end