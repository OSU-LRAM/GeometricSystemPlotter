function [scaled, min_coefficient_scaled] = arctan_scale_vector_fields(fields)
% Apply a scaled-arctan normalization to the vector fields that are the
% rows of field_cell. the second output uses the smallest scale out of the
% different fields, rather than the individual scale for each.


	% get the number of vector fields
	n_rows = size(fields,1);
	
	% preallocate output
	scaled = cell(size(fields));
	min_coefficient_scaled = scaled;
	magV = cell(size(fields,1),1);
	scalemag = zeros(size(fields,1),1);
	
	%%%%%%%%
	% Operate over each row
	for i = 1:n_rows
		
		%%%
		% Get the magnitude of each field
		
		% Concatenate the fields along a higher dimension
		n_dim = numel(size(fields{i,1}));
		
		magV{i} = sqrt(sum(cat(n_dim+1,fields{i,:}).^2,n_dim+1));
		
		% Get the scale factor
		scalemag(i) = .2/median(magV{i}(:)); %.2 is about the middle of the linear region of arctan
		if ~isfinite(scalemag(i))
			scalemag(i) = 1;
		end
		
		% Apply the arctan scaling
		scaled(i,:) = cellfun(@(x) atan_scaled(scalemag(i), x, magV{i}),fields(i,:),'UniformOutput',false);
		
	end
	
	%%
	% With all the scalemags found, scale using the smallest
	minscalemag = min(scalemag);
	for i = 1:n_rows
		
		min_coefficient_scaled(i,:) = cellfun(@(x) atan_scaled(minscalemag, x, magV{i}),fields(i,:),'UniformOutput',false);
		
	end



end


% Scaled arctangent function
function arrayout = atan_scaled(scalemag, arrayin, magV)

	arrayout = arrayin .* atan(scalemag*magV)./(scalemag*magV);
	
	arrayout( ~isfinite(arrayout) | isnan(arrayout) ) = 0;

end

