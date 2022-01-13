function [maxpoint, maxplane, Imax_x, Imax_y, Imax_z,interpstatecurvature] = CCF_maxpoint(CCF,grid)                    
    n_dim = numel(size(CCF{1}));

    % create an empty cell array with as many entries as
    % there are dimensions
    y=cell(1,n_dim);

    % populate this cell array with zeros
    y(:)={0};

    % Build up a (squared) norm of the constraint curvature
    % (using for now the simple in-coordinates norm, could
    % introduce the metric later)
    CCFnorm = zeros(size(CCF{1}));
    for idx_normbuild = 1:size(CCF,2)
        CCFnorm = CCFnorm+CCF{idx_normbuild}.^2;
    end


    
    
    % Get the maximum value, and the index of that value,
    % for the current CCF
    [~, Imax] = max(CCFnorm(:));
    [Imax_x,Imax_y,Imax_z] = ind2sub(size(CCF{1}),Imax);

    % Populate y with the location of that point
    for idx_maxpoint = 1:n_dim
        y{idx_maxpoint} = grid{idx_maxpoint}(Imax);
    end

    

%     % create a cell array with two fewer elements than
%     % there are dimensions
%     idxt=cell(1,n_dim-2);
% 
%     % populate this cell array with ones
%     idxt(1,:)={1};

    % For every dimension in the problem, create a copy of
    % the grid? Not sure why this is named
    % "interpstatecurvature"
    for j=1:1:n_dim
        interpstatecurvature{j}=grid{j,1};
    end

    % Get the value of the curvature at the point y
    for j=1:n_dim*(n_dim-1)/2
        curvature(:,j)=interpn(interpstatecurvature{:},CCF{j},y{:},'cubic');
    end

    % Create a skew-symmetric matrix describing the
    % orientation of the curvature form
    B=[0,curvature(1),curvature(2);
        -curvature(1),0,curvature(n_dim);
        -curvature(2),-curvature(n_dim),0];

    % Get the eigenvectors and values of the curvature
    % orientation
    [V,D]=eig(B);

    % Sort the eigenvalues and eigenvectors
    [d,ind] = sort(diag(D));
    Ds = D(ind,ind);
    Vs = V(:,ind);


    % Get the real component of the eigenvector with the
    % largest eigenvalue
    X=real(Vs(:,end));

    % Get a vector that is right-hand-positive orthogonal
    % to the preceding vector
    Y=(Vs(:,end)-X)/(sqrt(-1));

    % Normalize the X and Y vectors (to make up for the
    % fact that they lost length when mapped to real values
    Xnorm=X/(norm(X));
    Ynorm=Y/(norm(Y));
    
    maxpoint = y;
    maxplane = [Xnorm,Ynorm];
    
end