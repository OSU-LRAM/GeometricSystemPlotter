function solution = ode_multistart(odefun,argfun,limits,start_t,start_y,options)
% This function is a wrapper for Matlab's ODE functions that allows for the
% specification of several windows of integration, restarting the
% integrator at each point. The starting independent variable can also be
% specified inside the range, in which case the integrator will be run both
% forwards and backwards from this point. The arguments are:  
%
% odefun: handle to the ODE integrator to use
% argfun: handle to the function to be integrated
% limits: an ascending or descending *row* vector from which TSPAN will be drawn
% start_t: the starting indpendent variable value.
% start_y: the value of the output at start_t.
% options: an ODE options structure
%
% The output is a solution function that gives y(t) for any point inside
% the range, and NaN and a warning for points outside.

	% create an empty options structure if not provided
	if ~exist('options','var')
		options = [];
	end

	%%%%%%%%%%%
	% First, make sure that start_t is included as a limit of integration, and
	% ensure orderedness of the limits
	all_limits = unique([limits,start_t]);

	% Separate into limits for moving forward and backward from the midpoint
	pos_limits = all_limits(all_limits >= start_t);
	pos_limit_stack = [pos_limits(1:end-1)' pos_limits(2:end)'];

	neg_limits = all_limits(all_limits <= start_t);
	neg_limits = neg_limits(end:-1:1);
	neg_limit_stack = [neg_limits(1:end-1)' neg_limits(2:end)'];


	%%%%%%%%%%%%%
	% Second, run the integrator over each section

	% Positive direction
	pos_sol = cell(size(pos_limit_stack,1),1); %cell array to store output structures
	init_y = start_y; %Initial condition, will be updated at the end of each subrun

	for i = 1:length(pos_sol);

		% get the solution structure
		pos_sol{i} = odefun(argfun,pos_limit_stack(i,:),init_y,options);

		% get the final value to start the next step
		init_y = deval(pos_sol{i},(pos_limit_stack(i,2)));

	end

	% Negative direction
	neg_sol = cell(size(neg_limit_stack,1),1); %cell array to store output structures
	init_y = start_y; %Initial condition, will be updated at the end of each subrun

	for i = 1:length(neg_sol);

		% get the solution structure
		neg_sol{i} = odefun(argfun,neg_limit_stack(i,:),init_y,options);

		% get the final value to start the next step
		init_y = deval(neg_sol{i},(neg_limit_stack(i,2)));

	end

	%%%%%%%%%%%%%
	% Concatenate the solutions and limit stacks from the negative end to the
	% positive end

	all_sol = cat(1,neg_sol(end:-1:1),pos_sol);
	all_limit_stack = cat(1,neg_limit_stack(end:-1:1,end:-1:1),pos_limit_stack);

	%%%%%%%%%%
	% Make the output function

	solution = @(t) cell2mat(arrayfun(@(t) multideval(all_sol, all_limit_stack, t(:).'), t,'UniformOutput',false));
	
end

function y = multideval(all_sol, all_limit_stack, t)
% Evaluate a piecewise differential equation defined by solution structures
% in all_sol, over domains in all_limit_stack, at t

	% Establish which domain t is in
	within = ( all_limit_stack(:,1) <= t ) & ( all_limit_stack(:,2) > t );
	
	% Special case to make sure that upper bound of range is included
	within(end) = ( all_limit_stack(end,1) <= t ) & ( all_limit_stack(end,2) >= t );
	
	% Evaluate the corresponding solution structure
	switch sum(within)
		
		case 0 % t is out of bounds
			
			y = NaN(size(deval(all_sol{1},(all_limit_stack(1,1)))));
			
		case 1 % t has been uniquely placed
			
			y = deval(all_sol{within},(t));
			
		otherwise % t could not be uniquely placed
			
			error('t could not be uniquely placed within the time sequence')
			
	end

end