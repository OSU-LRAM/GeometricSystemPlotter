function M_alpha = mass_pull_back(M_full,A)

%     if ~isa(A,'cell')
%         A_calc = {A};
%         convertedA = 0;
%     else
%         A_calc = celltensorconvert(A);
%         convertedA = 1;
%     end
        
%     M_alpha = cell(size(A_calc));

%     for i = 1:numel(A_calc)
%         % Pull M_full back to be in terms of the shape variables
%         M_alpha{i} = [-A_calc{i}' eye(size(A_calc{1},2))]*M_full*[-A_calc{i}; eye(size(A_calc{1},2))];
%     end
%     
%     if convertedA
%         M_alpha = celltensorconvert(M_alpha);
%     end
    M_alpha = [-A' eye(size(A,2))]*M_full*[-A; eye(size(A,2))];
    
end