%Evaluate the connection over the fine grid for calculations and the coarse
%grid for vector display
function s = evaluate_connection(s)

	% If there is no connection denominator, make it all ones
	if ~isfield(s,'A_den')
		n_shape = nargin(s.A_num);
		shape_test_list = num2cell(ones(n_shape,1));
		s.A_den = @(a,varargin) repmat(ones(size(a)),size(s.A_num(shape_test_list{:})));
	end

    %list of zoom levels at which to evaluate connection vector fields and the
    %grid which should be used for them
    vector_field_list = {'display','vector';
                         'eval','eval'};
                     
    %loop over list, creating the vector fields
    for i = 1:size(vector_field_list,1);
        
        %generate one large array with all the vector information in it for
        %each of the numerator and denominator, and treat the set as two
        %cells in an array
        %note that the fields are the _negative_ of the connection
        s.vecfield.(vector_field_list{i,1}).content.A_num =... %numerator
            -s.A_num(s.grid.(vector_field_list{i,2}){:});
        s.vecfield.(vector_field_list{i,1}).content.A_den =... %denominator
            s.A_den(s.grid.(vector_field_list{i,2}){:});
        
        %mark what grid was used to create the field
        s.vecfield.(vector_field_list{i,1}).type = (vector_field_list{i,2});
        
    end

end