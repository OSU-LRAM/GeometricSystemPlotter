function output = shchf_diffdrive_cyclic(input_mode,pathnames)
	% Default argument
    if ~exist('input_mode','var')
		input_mode = 'initialize';
    end
    
    switch input_mode
		case 'name'
			output = 'Diffdrive Cyclic Motion';
		case 'dependency'
			output.dependency = {};
		case 'initialize'
            % path definition
            % equivalent to example motion, but cyclic (parallel parking)
            %side_lengths = pi/32:pi/64:pi/8;
            side_lengths = pi/24:pi/24:pi/6;
            p.phi_def = cell(1,4*length(side_lengths));
            ctr = 0;
            for i = 1:length(side_lengths)
                % phase
                for j = 0:3
                    ctr = ctr + 1;
                    cycle = {@(t) drive_forward(t, side_lengths(i)),...
                             @(t) turn_left(t, side_lengths(i)),...
                             @(t) drive_backward(t, side_lengths(i)),...
                             @(t) turn_right(t, side_lengths(i))};
                    p.phi_def{ctr} = circshift(cycle, j);
                    % enable area integration
                    p.cBVI_method{ctr}{1} = 'simple';
                end
            end
            
			%marker locations
			p.phi_marker = []; % No marker on this path (can put, e.g. endpoints of path if desired)
			
			%arrows to plot
			p.phi_arrows = 1;

			%time to run path
			p.time_def = [0 1]; % Duration of each segment

			%number of points in each path.
			p.phi_res = 20;

			% output the path properties
			output = p;
    end
end

function [alpha] = drive_forward(t,len)
	t = t(:)*len;
    start = [0 -len];
	alpha = start + [t t];
end

function [alpha] = turn_left(t,len)
	t = t(:)*len;
	alpha = drive_forward(1,len) + [-t t];
end

function [alpha] = drive_backward(t,len)
    t = t(:)*len;
    alpha = turn_left(1,len) + [-t -t];
end

function [alpha] = turn_right(t,len)
    t = t(:)*len;
    alpha = drive_backward(1,len) + [t -t];
end