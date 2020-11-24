function T = evaluate_tensor_over_grid(tensorfunction,grid,ignore_singular_warning)

% Check what type of output the tensorfunction produces in response to a
% grid input (block format in which each section is one component evaluated
% over the grid, block format in which each tile is the full tensor
% evaluated at a point, or gridded inputs result in function failure or
% ambiguity

%%%%%%%
% Allow user to specify that singular matrix warnings can be ignored
if ~exist('ignore_singular_warning','var')
    ignore_singular_warning = 0;
end


% Extract the grid range from the grid
grid_range = cellfun(@(g) [min(g(:)) max(g(:))],grid,'UniformOutput',false);
grid_range = cell2mat(grid_range(:)');

% Test the output of the function
tensorfunctiontype = test_function_type(tensorfunction,grid_range,ignore_singular_warning);



%%%%%%%%%%%%

% Either call the connection function directly for the shape
% inputs, or run a parallel loop to evaluate it at each of the
% points
switch tensorfunctiontype

    case 'block'

        % Evaluate the function at all points
        A_full = tensorfunction(grid{:});

        % Split the function output into pieces each the
        % size of the input
        T = mat2tiles(A_full, size(grid{1}));

    case 'woven'

        % Evaluate the function at all points
        A_full = tensorfunction(grid{:});

        % Split the function output into pieces each the
        % size of the output
        a_point = cellfun(@(ai) ai(1),grid,'UniformOutput',false);
        A_cell = mat2tiles(A_full, size(tensorfunction(a_point{:})));


        % Swap the inner and outer structure of the cell so
        % that outer layer is position in the tensor and
        % inner layer is position in the grid.
        T = celltensorconvert(A_cell);

    case 'single point'


        A_cell = cell(size(grid{1})); % Build a cell array to hold the function at each point
        parfor par_idx = 1:numel(grid{1})   % Loop over all elements of the grid
            
            if ignore_singular_warning
                warning('off','MATLAB:singularMatrix');
                warning('off','MATLAB:illConditionedMatrix');
            end

            % Extract the idx'th element of each grid, and put them
            % into a cell array
            a_point = cellfun(@(ai) ai(par_idx),grid,'UniformOutput',false);

            % Evaluate the function at the idx'th point
            A_cell{par_idx} = tensorfunction(a_point{:});
            
            warning('on','MATLAB:singularMatrix');
            warning('on','MATLAB:illConditionedMatrix');
            
        end

        % Swap the inner and outer structure of the cell so
        % that outer layer is position in the tensor and
        % inner layer is position in the grid.
        T = celltensorconvert(A_cell);

    case 'cell block'
        
        T = tensorfunction(grid{:});

    case 'cell woven'
        
        T = celltensorconvert(tensorfunction(grid{:}));
        
    otherwise

        error('Function type was specified as something unexpected')

end


end