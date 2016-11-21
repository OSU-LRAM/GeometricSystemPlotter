% Calculate the tension in each spring. A stretched spring has positive
% tension. Using a logarithmic force rule for geometric scale invariance,
% with a spring constant of 1. (the
% ratio of the spring constant to the first-order drag coefficient gives
% the timescale, but we don't really care about the timescale)
function T = get_spring_tensions(current_length,neutral_length,masked_springs)

    % Force in a logarithmic spring
	T = log((current_length./neutral_length));
    
    % Zero out the force on any spring outside the masked region
    T = T.*masked_springs;
	
end
