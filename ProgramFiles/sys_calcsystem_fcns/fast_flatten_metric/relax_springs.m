%Integrate to relax the springs
function [final_x,final_y, sol] = relax_springs(x,y,springs,neutral_lengths,tol)
	%tol is allowable ratio of movement between iterations and from first to last
	
	%Set up the ODE for moving the nodes
	IC_start = [x(:);y(:)];
	IC = IC_start;
	
	% Iterate doubling integration lengths until the difference between the
	% most recent run and the
		
    time_span = 1; %Check ODE every this often
	    
	% Fast solution to prime automated convergence test
	sol_last = ode45(@(t,state) relaxation_helper_function(t,state,neutral_lengths,springs),[0 time_span/10000],IC);
    
	
    last_x = x(:);
    last_y = y(:);
    current_x = x(:);
    current_y = y(:);
    last_dist2_moved = inf;
    
    counter = 1;
    slowing = 0; % Use this to keep track of inflection on distance moved
	while 1
		
        time_limit = counter*time_span;
		% Integrate the spring relaxation
		sol_test = odextend(sol_last, @(t,state) relaxation_helper_function(t,state,neutral_lengths,springs),time_limit);
		
        % Get the current positions of the nodes
        current_state = deval(sol_test,counter*time_span);
        last_state = deval(sol_test,(counter-1)*time_span);
%         current_x(:) = current_state(;
%         current_y(:) = sol_test.y((end/2 + 1):end,end);
       

        % Get the squared distance the points have moved
        current_dist2_moved = sum((current_state-last_state).^2);%sum( (current_x(:)-last_x(:)).^2 + (current_y(:)-last_y(:)).^2 );
        
        % Get current lengths
        current_lengths = get_spring_lengths_and_azimuths(springs,current_x,current_y);
        
        % Get current length ratios
        current_length_ratios = current_lengths(:)./neutral_lengths(:);
        
        % Check termination condition: Run until springs have all settled
        % to within tol of their neutral lengths, or until they have
        % started moving faster after having slowed down
        
        if ~slowing && counter ~=1
            if current_dist2_moved < .5*last_dist2_moved;
                slowing = 1;
            end
        end
        
        if (current_dist2_moved > last_dist2_moved) && slowing
            sol = sol_last;
            break
        elseif max(abs(log(current_length_ratios))) < log(1+tol)
            sol = sol_test;
            break
        else
            last_dist2_moved = current_dist2_moved;
            sol_last = sol_test;
            last_x = current_x;
            last_y = current_y;
        end
        
        
		
		%If the safety counter has run out, stop integrating
		if counter > 100;
			
            sol = sol_test;
			break			
			
		else
						
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