%Reshape the evaluated metric to format compatible with plotting over
%base space
function s = reshape_metric(s)

	%get the list of metric field zoom levels
	zoom_list = fieldnames(s.metricfield);


	%loop over all zoom levels of the field being evaluated
	for m = 1:length(zoom_list)

		%loop over all contributing parts of the connection
		contributors = fieldnames(s.metricfield.(zoom_list{m}).content);

		for n = 1:length(contributors)

			%move the flat version out of the way into local storage
			metricfield_orig = ...
				s.metricfield.(zoom_list{m}).content.(contributors{n});


			%get the target size of the field being evaluated
			target_size = size(s.metricfield.(zoom_list{m}).content.(contributors{n}){1});

			%clear out room for reshaped version
			s.metricfield.(zoom_list{m}).content.(contributors{n}) = cell(target_size);



			%loop over each row of the conection
			for i = 1:prod(target_size)

				s.metricfield.(zoom_list{m}).content.(contributors{n}){i} = ...
					cellfun(@(a) a(i),metricfield_orig);
				

			end

		end



	end


end