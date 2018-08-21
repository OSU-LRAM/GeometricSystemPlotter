%function to divide numerator/denominator, and deal with singularities 

function [quotient, flag_array_out] = smart_divider(grid,nums,dens)

	%%%%%%%%%%%%%
	%Main dividing code

	%divide numerator by denominator
	quotient = nums./dens;
	flag_array_out = zeros(size(nums)); %to hold singularities

	%small number to replace zeros in denominator with and big number to
	%set singularities to if necessary
	smallnum = 1e-8;
	%bignum = 1/smallnum;



	%find all entries with zeros in the denominator
	zero_indices = ( (dens > -smallnum) & (dens < smallnum) );

	% Identify the singularities

	flag_array_out(zero_indices) = 1;

	% Handle any singularities identified
	if any(zero_indices)
		
		quotient(zero_indices) = 0;

% 		% Get the number of dimensions
% 		n_dim = numel(size(nums));
% 
% 		% Turn the grid into a set of columns
% 		gridcolumns = grid_to_columns(grid);
% 
% 		% Use TriScatteredInterp to fill in singularities if dimensionality
% 		% is small enough
% 		if n_dim <=3
% 
% 			% Get the interpolating function
% 			f_fix = TriScatteredInterp(gridcolumns(~zero_indices,:),quotient(~zero_indices),'natural');
% 
% 			% Fill in the singularities
% 			quotient(zero_indices) = f_fix(gridcolumns(zero_indices,:));
% 
% 			% Otherwise, use griddata
% 		else
% 
% 			quotient(zero_indices) = griddatan(gridcolumns(~zero_indices,:),quotient(~zero_indices));
% 
% 		end

	end
	
end

