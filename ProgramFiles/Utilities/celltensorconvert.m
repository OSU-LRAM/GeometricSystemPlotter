function B = celltensorconvert(A)
% Takes an NxMx... cell array of YxZx... arrays, and returns an YxZx... cell array of
% NxMx.. arrays


% Create a cell containing a matrix whose dimensions are taken from A
new_cell = {zeros(size(A))};

% Replicate this cell to make a cell array whose dimensions are taken from
% the contents of (the first value of) A
B = repmat(new_cell,size(A{1}));

for i = 1:numel(A)
	
	for j = 1:numel(B)
		
		B{j}(i) = A{i}(j);
		
	end
	
end
