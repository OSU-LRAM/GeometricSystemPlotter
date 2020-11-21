function [s,T,T2] = evaluate_tensor_over_grid(s,tensorfunction,grid,ignore_singular_warning,A_eval,A_grid)

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

if exist('A_eval','var') && ~isempty(A_eval)
    super_special_condition = 1;
else
    super_special_condition = 0;
end

%Initialize second tensor for return in case we don't set it later
T2 = [];

% Extract the grid range from the grid
grid_range = cellfun(@(g) [min(g(:)) max(g(:))],grid,'UniformOutput',false);
grid_range = cell2mat(grid_range(:)');

% Test the output of the function
if super_special_condition
    tensorfunctiontype = test_function_type(s,tensorfunction,grid_range,ignore_singular_warning,A_eval,A_grid);
else
    tensorfunctiontype = test_function_type(s,tensorfunction,grid_range,ignore_singular_warning);
end


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
        if isequal(tensorfunction,s.A_num)
            Ma_cell = cell(size(grid{1}));
            M_full_cell = cell(size(grid{1}));
            J_full_cell = cell(size(grid{1}));
            local_inertias_cell = cell(size(grid{1}));
        end
        
        for par_idx = 1:numel(grid{1})   % Loop over all elements of the grid
            
            if ignore_singular_warning
                warning('off','MATLAB:singularMatrix');
                warning('off','MATLAB:illConditionedMatrix');
            end

            % Extract the idx'th element of each grid, and put them
            % into a cell array
            a_point = cellfun(@(ai) ai(par_idx),grid,'UniformOutput',false);
            a_vec = cell2mat(a_point');

            % Evaluate the function at the idx'th point
            if super_special_condition
                try
                    J_full = s.evaluated_vals.J_full{par_idx};
                    local_inertias = s.evaluated_vals.local_inertias{par_idx};
                    M_full = s.evaluated_vals.M_full{par_idx};
                    T = tensorfunction(M_full,J_full,local_inertias,cell2mat(a_point));
                catch
                    T = tensorfunction(a_point{:},A_eval,A_grid);
                end
            else
                if isequal(tensorfunction,s.A_num)
                    [A, M_a,J_full, local_inertias,M_full] = tensorfunction(a_point{:});
                    T = A;
                else
                    T = tensorfunction(a_point{:});
                end
            end
            A_cell{par_idx} = T;
            if isequal(tensorfunction,s.A_num)
                Ma_cell{par_idx} = M_a;
                M_full_cell{par_idx} = M_full;
                J_full_cell{par_idx} = J_full;
                local_inertias_cell{par_idx} = local_inertias;
            end
            
            warning('on','MATLAB:singularMatrix');
            warning('on','MATLAB:illConditionedMatrix');
            
        end
        
        if isequal(tensorfunction,s.A_num)
            T2 = celltensorconvert(Ma_cell);
            s.evaluated_vals.M_full = M_full_cell;
            s.evaluated_vals.J_full = J_full_cell;
            s.evaluated_vals.local_inertias = local_inertias_cell;
        end

        if super_special_condition
            T_temp = cell_layer_swap(A_cell);
            T = cell(size(T_temp));
            for i = 1:length(T_temp)
                T{i} = celltensorconvert(T_temp{i});
            end
        else
            % Swap the inner and outer structure of the cell so
            % that outer layer is position in the tensor and
            % inner layer is position in the grid.
            T = celltensorconvert(A_cell);
        end

    case 'cell block'
        
        T = tensorfunction(grid{:});

    case 'cell woven'
        
        T = celltensorconvert(tensorfunction(grid{:}));
        
    otherwise

        error('Function type was specified as something unexpected')

end


end