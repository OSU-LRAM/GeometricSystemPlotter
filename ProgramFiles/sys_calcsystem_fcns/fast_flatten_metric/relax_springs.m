%Integrate to relax the springs
function [final_x,final_y, sol] = relax_springs(x,y,springs,neutral_lengths,tol)
	%tol is allowable ratio of movement between iterations and from first to last
	
	%Set up the ODE for moving the nodes
	IC_start = [x(:);y(:)];
	IC = IC_start;
	
	% Iterate doubling integration lengths until the difference between the
	% most recent run and the
		
	time_limit = 10; % Base time limit
	
	% Fast solution to prime automated convergence test
	sol = ode45(@(t,state) relaxation_helper_function(t,state,neutral_lengths,springs),[0 time_limit/10],IC);
	
    current_x = x;
    current_y = y;
    
    counter = 1;
	while 1
		
		% Integrate the spring relaxation
		sol = odextend(sol, @(t,state) relaxation_helper_function(t,state,neutral_lengths,springs),time_limit);
		
% 		% Get the change since the last run, and the change since the start
% 		output_change = sqrt(sum((sol.y(:,end)-ytest).^2));
% 		total_change = sqrt(sum((sol.y(:,end)-IC_start).^2));
        
        % Get the current lengths of the springs
        current_x(:) = sol.y(1:(end/2),end);
        current_y(:) = sol.y((end/2 + 1):end,end);
        current_lengths = get_spring_lengths_and_azimuths(springs,current_x,current_y);
        
        % Get the current velocities of the points
        current_velocities = relaxation_helper_function(0,sol.y(:,end),neutral_lengths,springs);
		
		%If the system has converged or safety counter has run out, stop integrating
		if (max(current_velocities) <= tol*min(current_lengths)) || counter > 6;
			
			break
			
		else
			
			%double the time limit
			time_limit = 2*time_limit;
			
			counter = counter + 1;
			
			
		end
	
		
	
	end

	
	% Extract the final x and y positions
	[final_x, final_y] = deal(zeros(size(x)));
	final_x(:) = sol.y(1:(end/2),end);
	final_y(:) = sol.y((1+end/2):end,end);
	
end

function dstate = relaxation_helper_function(~,state,neutral_lengths,springs)

	% Calculate the spring lengths and headings
	[lengths, az] = get_spring_lengths_and_azimuths(springs,state(1:(length(state)/2)),state(((length(state)/2)+1):end));
	
	% Find the tensions in the springs
    masked_springs = springs(:,4);
	T = get_spring_tensions(lengths,neutral_lengths,masked_springs);
	
	% Calculate the forces on the nodes
	[Fx,Fy] = get_node_forces(springs,T,az,[length(state)/2, 1]);
	
	% The forces specify the state velocity
	dstate = [Fx;Fy];

end