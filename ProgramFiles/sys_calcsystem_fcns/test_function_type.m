% Test whether the a function can be called on a grid of points, and
% if so, whether it produces block- or woven-style output
function tensorfunctiontype = test_function_type(tensorfunction,grid_range,ignore_singular_warning)

%%%%%%%
% Allow user to specify that singular matrix warnings can be ignored
if ~exist('ignore_singular_warning','var')
    ignore_singular_warning = 0;
end
if ignore_singular_warning
    warning('off','MATLAB:singularMatrix');
    warning('off','MATLAB:illConditionedMatrix');
end



% Extract the centerpoint of the range of the grid
range_start = grid_range(1:2:end);
range_end = grid_range(2:2:end);

range_mid = num2cell((range_start+range_end)/2);


% Try to evaluate the local connection numerator at this midpoint, but
% supply the point twice as if it were two grid points
double_midpoint = mat2tiles([range_mid{:} range_mid{:}],size(range_mid));
vector_input = 1;
A_test_single = tensorfunction(range_mid{:}); % This is to make the system throw an error if tensorfunction is bad, because the try line on double_midpoint will treat an error as indicating a non-vector function
try

    A_test = tensorfunction(double_midpoint{:});

catch
    
    % If the system fails the test, then the (fieldname) function can only take a
    % single point at a time
    
    vector_input = 0;
    
    
end

% If the system failed the vector input test, then mark it as a
% single-point evaluation
if  ~vector_input    

    tensorfunctiontype = 'single point';
    

elseif iscell(A_test)
    
    if all(size(A_test{1}) == size(double_midpoint{1}))
        
        tensorfunctiontype = 'cell block';
        
    elseif all(size(A_test) == size(double_midpoint{1}))
        
        tensorfunctiontype = 'cell woven';
        
    else
        
        error([func2str(tensorfunction) ' returns a cell array that doesn''t match a recognized tensor format'])
        
    end
        
% If the output is not a cell array, and A_test does not have twice as many
% entries as A_test_single, then this function doesn't handle vector inputs
% properly, and should be marked as single
elseif numel(A_test) ~= 2*numel(A_test_single)
    
    tensorfunctiontype = 'single point';
    
% If the function appears able to take vector input and return non-cell output, test to see if the
% output is block or woven
else
    
    % Evaluate the connection numerator at the midpoint, or at a random other
    % point if the cell/woven nature is ambiguous.

    counter = 0;
    ambiguous = 1;
    while 1

        % Evaluate connection numerator at midpoint
        A_mid = tensorfunction(range_mid{:});


        % Tile this out into two copies side by side
        A_mid_woven = [A_mid A_mid];


        % Make a matrix in which each block is one entry of (fieldname) tiled out
        A_mid_block_cell = arrayfun(@(x) x*[1 1],A_mid,'UniformOutput',false);
        A_mid_block = cell2mat(A_mid_block_cell);


        % Check to be sure that A_mid_block and A_mid_woven are different
        if all(A_mid_block(:) == A_mid_woven(:))

            % Midpoint test does not reveal any difference between block and
            % woven format

            % Generate a random point within the domain
            range_mid = num2cell(range_start + (0.9 * (range_end-range_start).*rand(size(range_end))));

            % Add one to the test counter
            counter = counter +1;

        else

            % Midpoint test has enough information to distinguish between woven
            % and block forms, so break out of loop

            ambiguous = 0;
            break;

        end

        % Stop trying more points after failing 10 tests.
        if counter == 10

            break

        end

    end

    % If the function's output is still ambiguous, treat it as a
    % single-point function. Otherwise, identify whether it is block or
    % woven.
    if ambiguous
  
        tensorfunctiontype = 'single point';
    
    else
    
        % Double-evaluate the connection at the midpoint or random test point
        % that was found above

        double_midpoint = mat2tiles([[range_mid{:}]' [range_mid{:}]'],size(range_mid));
        A_test = tensorfunction(double_midpoint{:});

        %%%%%%%%
        %%%%%%%%
        % Compare A_test with A_mid_block and A_mid_woven
        if all(A_test(:) == A_mid_woven(:))

            tensorfunctiontype = 'woven';

        elseif all(A_test(:) == A_mid_block(:))

            tensorfunctiontype = 'block';
            
        else
            
            tensorfunctiontype = 'single point';
            warning('check  the test_function_type program to make sure it is handling this function correctly')
            
        end

    end

end
    
warning('on','MATLAB:singularMatrix');
warning('on','MATLAB:illConditionedMatrix');


    
end
