% calc_ccf.m
% computes constraint curvature over a grid
% inputs:
    % grid_points: cell array of point arrays over which A is defined
    % A: local connection tensor (any position, shape space)
    % mat_fun:
        % anonymous function, mapping configuration vector representation
        % to corresponding matrix representation
    % vec_fun:
        % inverse of mat_fun; should return column vector
% outputs:
    % DA: total lie bracket (CCF) vector
    % dA: exterior derivative of A
        % represented as (n_dim)x(n_2_forms) cell array
        % involves choice of basis order - currently (earlier)^(later)
        % todo: order specification requires:
            % different indexing method
            % structure containing signs
    % lb: local Lie bracket of A columns

function [DA, dA, lb] = calc_ccf(grid_points, A, mat_fun, vec_fun)
    %% compute dA 2-form
    % preallocate dA: n_dim(xyz) x shvar_dim choose 2 (num. of 2-forms)
    dA = cell(size(A,1), nchoosek(size(A,2), 2));
    % iter. through ea. dim. (xyz) in A
    for i = 1:size(A,1)
        % compute gradient of ea. shape w.r.t. ea. other
        % repr. as J (rows are shvars, cols are der. wrt. other shvar)
        partials = cell(size(A,2));
        order = [2 1 3:size(A,2)]; %ordering for gradient
        grad_Aij = cell(1, size(A,2)); %temp. storage for ea. gradient
        for j = 1:size(A,2)
            % get gradient (partial of this shape wrt ea. other, everywhere)
            [grad_Aij{order}] = gradient(A{i,j},grid_points{order});
            partials(j,:) = grad_Aij;
        end
        % iter. thru J (all rows, col < row):
        entry = 1;
%         for r = 1:size(A,2)
%             for c = 1:r
        for c = 1:size(A,2)
            for r = c:size(A,2)
                % combine like bases (respecting wedge product rules)
                % if doing a different choice of basis, change:
                    % which entry (r,c) pairs are saved in
                    % sign choices for TR and BL off-diagonals
                if r == c
                    % skip zero wedge products
                    continue;
                end
                % two partials for this basis
                dA{i,entry} = partials{r,c} - partials{c,r};
                
                entry = entry + 1;
            end
        end
    end
    
    %% compute local Lie bracket
    % convert A; working pointwise
    A_pts = celltensorconvert(A);
    % compute Lie bracket at ea. point (parfor -> speed)
    lb_pts = cell(size(A_pts));
    parfor i = 1:numel(A_pts)
        % suppressed warning about function calls in parfor; unavoidable
        %#ok<*PFBNS>
        % preallocate space
        lb_pts{i} = zeros(size(dA));
        % iterate thru. combinations of bases
        entry = 1;
%         for aj = 1:size(A,2)
%             for ak = 1:aj
        for ak = 1:size(A,2)
            for aj = ak:size(A,2)
                % skip zero wedge products
                if aj == ak
                    continue;
                end
                % get local Lie bracket of shapes
                % done according to configuration space
                X = mat_fun(A_pts{i}(:,aj));
                Y = mat_fun(A_pts{i}(:,ak));
                lb_pts{i}(:,entry) = vec_fun(X*Y-Y*X);
                entry = entry + 1;
            end
        end
    end
    lb = celltensorconvert(lb_pts);
    
    %% combine for DA
    DA = cell(size(dA));
    for i = 1:numel(dA)
        DA{i} = dA{i} - lb{i}; %sign choice for consistency
    end
 
end