function [gradE, E, resid] = fasthelmholtz(grid,V)

% Recreate the original spacing vectors defining the grid
gridSpacing = extractGridSpacing(grid);

% % Convert from ndgrid to meshgrid
% dimcount = numel(grid);     % Count how many dimensions field is defined over
% dimlist = [2 1 3:dimcount]; % Permutation mapping from ndgrid to meshgrid
% grid = permute(grid,dimlist); % Perform permutation on grid component ordering
% V = permute(V,dimlist);       % Perform permutation on vector component ordering 
% gridSpacing = permute(gridSpacing,dimlist); % Perform permutation on grid spacing component ordering
% %grid = cellfun(@(g) permute(g,dimlist),grid,'UniformOutput',false); % Perform permutation on grid value ordering
% %V = cellfun(@(g) permute(g,dimlist),V,'UniformOutput',false);       % Perform permutation on vector value ordering


grid = {grid{2}; grid{1}; grid{3:end}};
V = {V{2}; V{1}; V{3:end}};


%%%%
% Take the gradient of each vector component

% Make a cell matrix the same size as the vector cell matrix
gradV = cell(size(V));

% Loop over each element of this cell (one cell per component of the
% vector)
for idx = 1:numel(grid)
    
    % Each vector component has a derivative along each vector component
    gradV{idx} = cell(size(grid));
    
    % Calculate the gradient
    [gradV{idx}{:}] = gradient(V{idx},gridSpacing{:});
end


%%%%%
% Convert the gradients into gradients calculated at each location

% One matrix at each grid point, with as many rows/columns as there are
% dimensions
gradV_local = repmat({zeros(numel(grid),numel(grid))},size(grid{1}));

for idx1 = 1:numel(grid)              %Loop over all vector components
    for idx2 = 1:numel(grid)          %Loop over all derivative directions
        for idx3 = 1:numel(grid{1})   %Loop over all grid points
            
            gradV_local{idx3}(idx1,idx2) = gradV{idx1}{idx2}(idx3);
            
        end        
    end
end


%%%%%%%%%%
% Remove the skew-symmetry in the gradient, remainder is derivative of
% conservative function

gradV_conservative = gradV;

for idx1 = 1:numel(grid)              %Loop over all vector components
    for idx2 = 1:numel(grid)          %Loop over all derivative directions
        for idx3 = 1:numel(grid{1})   %Loop over all grid points
            
            gradV_conservative{idx1}{idx2}(idx3) = ...
                (gradV{idx1}{idx2}(idx3) + ...
                    gradV{idx2}{idx1}(idx3))/2;
                
        end
    end
end


%%%%%%%
% Integrate conservative component
V_conservative = V;
for idx = 1:numel(grid)   % Loop over all vector components
    V_conservative{idx} = intgrad2(gradV_conservative{idx}{:},gridSpacing{:});
end

%%%%%%%
% Subtract conservative component from original vector field to get the
% residual non-conservative component
V_nonconservative = V;
for idx = 1:numel(grid)   % Loop over all vector components
    V_nonconservative{idx} = V{idx} - V_conservative{idx};
end

%%%%%%%
% Calculate the mean of the nonconservative field, so that we can extract
% this harmonic component
% residual non-conservative component
V_nonconservative_mean = V;
for idx = 1:numel(grid)   % Loop over all vector components
    V_nonconservative_mean{idx} = mean(V_nonconservative{idx}(:));
end

%%%%%%%%
% Add the mean of the non-conservative field back into the conservative
% field
V_conservative_fixed = V;
V_nonconservative_fixed = V;
for idx = 1:numel(grid)   % Loop over all vector components
    V_conservative_fixed{idx} = V_conservative{idx}+V_nonconservative_mean{idx};
    V_nonconservative_fixed{idx} = V_nonconservative{idx}-V_nonconservative_mean{idx};
end


% 
% %%%%
% % Split out the conservative component of the vector gradient
% 
% % Predefine a cell array to hold conservative components
% gradV_local_conservative = gradV_local;
% 
% for idx = 1:numel(gradV_local)  % Loop over all grid points
%     
%    % Conservative component is average of gradient and its transpose
%    gradV_local_conservative{idx} = (gradV_local{idx} + transpose(gradV_local{idx}))/2;
%     
% end


%%%%%
% Make the 


% % Convert from meshgrid to ndgrid
V_conservative_fixed = {V_conservative_fixed{2}; V_conservative_fixed{1}; V_conservative_fixed{3:end}};
% V_conservative_fixed = cellfun(@(g) permute(g,dimlist),V_conservative_fixed,'UniformOutput',false);       % Perform permutation on vector value ordering

V_nonconservative_fixed = {V_nonconservative_fixed{2}; V_nonconservative_fixed{1}; V_nonconservative_fixed{3:end}};
% V_nonconservative_fixed = cellfun(@(g) permute(g,dimlist),V_nonconservative_fixed,'UniformOutput',false);       % Perform permutation on vector value ordering

E = [];

gradE = V_conservative_fixed;
resid = V_nonconservative_fixed;


end