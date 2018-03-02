function robot = draw_swimmer(robot,suppress_line,line_at_origin)

robot.links = draw_set(robot.links);

if ~exist('suppress_line','var') || (suppress_line == 0)

	%Put the center-dot at the ref position
	%set(robot.center,'XData',robot.ref_pos(1),'YData',robot.ref_pos(2),'Zdata',3)

	dashed_line = [0:.2:2; (0:.2:2)+.12; NaN(size(0:.2:2))]/6;

	dashed_line_x = dashed_line(:);
	dashed_line_y = zeros(size(dashed_line_x));

	dashed_line_xy = ([cos(robot.ref_pos(3)) -sin(robot.ref_pos(3)); sin(robot.ref_pos(3)) cos(robot.ref_pos(3))] * [dashed_line_x'; dashed_line_y'])';

	dashed_line_xy = [-dashed_line_xy(end:-1:1,:); dashed_line_xy([1,end-1],:)];

	
	if ~exist('line_at_origin','var') || (line_at_origin == 0)
		linecx = robot.ref_pos(1);
		linecy = robot.ref_pos(2);
	else
		linecx = 0;
		linecy = 0;
	end
	
	set(robot.orientation,'XData',dashed_line_xy(:,1)+linecx,'YData',dashed_line_xy(:,2)+linecy,'ZData',2.5*ones(size(dashed_line_xy(:,1))));

end

end