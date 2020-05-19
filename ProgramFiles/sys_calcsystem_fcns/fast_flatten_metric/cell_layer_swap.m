function B = cell_layer_swap(A)
% Takes an NxM cell array of YxZ cell arrays, and returns an YxZ cell array of
% NxM cell arrays


% Create a cell containing a cell whose dimensions are taken from A
new_cell = {cell(size(A))};

% Replicate this cell to make a cell array whose dimensions are taken from
% the contents of (the first value of) A
B = repmat(new_cell,size(A{1}));

for i = 1:numel(A)
	
	for j = 1:numel(B)
		
		B{j}{i} = A{i}{j};
		
	end
	
end
