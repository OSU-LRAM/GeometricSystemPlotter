%If connection is defined with separate numerator and denominator, merge
%them. In any case, make Avec the *negative* of the local connection.

function s = merge_connection(s)

    %get the list of connection vector field zoom levels
    zoom_list = fieldnames(s.vecfield);
    
%     %Get the size of the local connection matrix
%     A_size = size(s.vecfield.(zoom_list{1}).content.A_num);
    
    %loop over all zoom levels of the field being evacuated
	for m = 1:length(zoom_list)
        
% 		% Allocate a cell array to hold the local connection
% 		s.vecfield.(zoom_list{m}).content.Avec = cell(A_size);

        % If no denominator, just copy the negative of A into Avec
        if ~isfield(s,'A_den')
            
            s.vecfield.(zoom_list{m}).content.Avec = ...
                cellfun(@(A) -A, s.vecfield.(zoom_list{m}).content.A_num,'UniformOutput',false);
            
        else % Divide in the denominator, but look for and handle singularities while doing this
            
        % The negative sign in the local connection is applied here.
        
		[s.vecfield.(zoom_list{m}).content.Avec ...
			,s.vecfield.(zoom_list{m}).singularities]...
			= cellfun(...
            @(num,den) smart_divider(s.grid.(s.vecfield.(zoom_list{m}).type),-num,den)... %note negative sign in num argument
			,s.vecfield.(zoom_list{m}).content.A_num...
			,s.vecfield.(zoom_list{m}).content.A_den...
            , 'UniformOutput',false);
		
	end

end