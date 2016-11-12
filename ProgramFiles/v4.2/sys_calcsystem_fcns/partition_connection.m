%Reshape the evaluated connection to format compatible with plotting over
%base space
function s = partition_connection(s)

	%get the list of connection vector field zoom levels
	zoom_list = fieldnames(s.vecfield);


	%loop over all zoom levels of the field being evaluated
	for m = 1:length(zoom_list)

		%loop over all contributing parts of the connection
		contributors = fieldnames(s.vecfield.(zoom_list{m}).content);

		for n = 1:length(contributors)

			%move the flat version out of the way into local storage
			vecfield_orig = ...
				s.vecfield.(zoom_list{m}).content.(contributors{n});

			% Get the size of the raw original field
			raw_size = size(vecfield_orig);

			% Build a cell array of 1s, for all dimensions of the raw
			% vector field
			input_dims = repmat({1},[numel(raw_size), 1]);

			% Build a cell array of colons, for all but the first two
			% dimensions of the raw vector field
			extra_dims = repmat({':'},[numel(raw_size)-2, 1]);

			%get the target size of the field being evaluated
			target_size = size(s.(contributors{n})(input_dims{:}));

			%clear out room for reshaped version
			s.vecfield.(zoom_list{m}).content.(contributors{n}) = cell(target_size);

			% Get the size of the chunks to break up along the first two
			% dimensions
			chunk_size = raw_size(1:2)./target_size;



			%loop over each row of the conection
			for i = 1:target_size(1)

				%loop over each column of the connection
				for j = 1:target_size(2)

					%get the rows and columns corresponding to the target patch
					rows = (i-1)*chunk_size(1) + (1:chunk_size(1));

					columns = (j-1)*chunk_size(2) +(1:chunk_size(2));

					%reshape the vector field
					s.vecfield.(zoom_list{m}).content.(contributors{n}){i,j} =...
						vecfield_orig(rows,columns,extra_dims{:});

				end

			end

		end



	end


end