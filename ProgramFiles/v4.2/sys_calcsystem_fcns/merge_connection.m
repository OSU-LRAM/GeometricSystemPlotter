%Merge components of evaluated connection
function s = merge_connection(s)

    %get the list of connection vector field zoom levels
    zoom_list = fieldnames(s.vecfield);
    
    %Get the size of the local connection matrix
    A_size = size(s.vecfield.(zoom_list{1}).content.A_num);
    
    %loop over all zoom levels of the field being evacuated
    for m = 1:length(zoom_list)
        
		% Allocate a cell array to hold the local connection
		s.vecfield.(zoom_list{m}).content.Avec = cell(A_size);
		
        %first, apply the smart divider to
        %get the base connection vector field Avec and the location of any
        %singularities
		[s.vecfield.(zoom_list{m}).content.Avec ...
			,s.vecfield.(zoom_list{m}).singularities]...
			= cellfun(@(num,den) smart_divider(s.grid.(s.vecfield.(zoom_list{m}).type),num,den)...
			,s.vecfield.(zoom_list{m}).content.A_num...
			,s.vecfield.(zoom_list{m}).content.A_den, 'UniformOutput',false);
		
	end

end