function vis = visibility(xp,yp)
% Build the visibility graph of points in the polygon described by xp and
% yp

	% Temporarily remove the endpoint if it is a duplication of the start
	% point
	if (xp(1) == xp(end)) && (yp(1) = yp(end))
		
		xp = xp(1:end-1);
		yp = yp(1:end-1);
		
		replace_end_flag = 1;
		
	else
		
		replace_end_flag = 0;
		
	end
	
	
	% Create the arrays of point-difference vectors, in which the i,jth
	% element is the relevant component of the vector from point i to point
	% j
	xdiff = cell2mat(arrayfun(@(p) xp(:)'-p,xp(:),'UniformOutput',false));
	ydiff = cell2mat(arrayfun(@(p) yp(:)'-p,yp(:),'UniformOutput',false));
	
	% Get the azimuth of each vector
	
	
	
	
end