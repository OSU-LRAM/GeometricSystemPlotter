function y = concatenate_functions(all_functions, all_limits, t)
% Evaluate a piecewise differential equation defined by solution structures
% in all_sol, over domains in all_limit_stack, at t

	% Get a permutation order to return the output with t in the right
	% direction
	if size(t,2) > size(t,1)
		output_order = [2 1];
	else
		output_order = [1 2];
	end

	y = arrayfun(@(T) permute(concatenate_helper(all_functions,all_limits,T),output_order)...
		,t,'UniformOutput',false);

	y = cell2mat(y);

end


function y = concatenate_helper(all_functions, all_limits, t)

	% Establish which domain t is in
	within = ( all_limits(:,1) <= t ) & ( all_limits(:,2) > t );
	
	% Special case to make sure that upper bound of range is included
	within(end) = ( all_limits(end,1) <= t ) & ( all_limits(end,2) >= t );
	
	% Evaluate the corresponding solution structure
	switch sum(within)
		
		case 0 % t is out of bounds
			
			y = NaN(size(all_functions{1}(all_limit_stack(1,1))));
			
		case 1 % t has been uniquely placed
			
			y = all_functions{within}(t);
			
		otherwise % t could not be uniquely placed
			
			error('t could not be uniquely placed within the time sequence')
			
	end

end