% Calculate the net x and y forces on each node
function [Fx, Fy] = get_node_forces(springs,T,start_az,gridsize)

	[Fx, Fy] = deal(zeros(gridsize));

	% Get the x and y forces exerted by springs
	spring_Fx = T.*cos(start_az);
	spring_Fy = T.*sin(start_az);
	
	% Apply the spring forces to the nodes
	for i = 1:size(springs,1)
		
		% Force applied to the base node
		Fx(springs(i,1)) = Fx(springs(i,1)) - spring_Fx(i);
		Fy(springs(i,1)) = Fy(springs(i,1)) - spring_Fy(i);
		
		% Force applied to the end node
		Fx(springs(i,2)) = Fx(springs(i,2)) + spring_Fx(i);
		Fy(springs(i,2)) = Fy(springs(i,2)) + spring_Fy(i);
		
	end

end
