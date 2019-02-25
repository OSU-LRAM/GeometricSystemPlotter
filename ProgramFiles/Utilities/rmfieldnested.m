function [s_out] = rmfieldnested(s,fieldarray)
%rmfieldnested: This function will return the structure s with the field at
%the end of the fieldarray removed. It will also remove any parent fields
%that are empty after doing this
% Inputs:
% s: struct to be removed from
% fieldarray: the heirarchy of fields to the one that should be deleted

s_out = s;

%%%%%%%%%
%First check to make sure that the field exists
if isfieldnested(s,fieldarray)
    if size(fieldarray,2) > 1
        %Loop through to remove that field
        currlevel{size(fieldarray,2)} = struct;
        currlevel{1} = s_out;
        for k = 1:size(fieldarray,2)-1
            currlevel{k+1} = currlevel{k}.(fieldarray{k});
        end
        
        %Remove the final element
        currlevel{end} = rmfield(currlevel{end},fieldarray{end});
        
        %Repopulate back into the original array
        for k = size(fieldarray,2)-1:-1:1
            currlevel{k}.(fieldarray{k}) = currlevel{k+1};
            %test = structfun(@isempty,currlevel{k}.(fieldarray{k}));
            if ~numel(fieldnames(currlevel{k}.(fieldarray{k})))
                currlevel{k} = rmfield(currlevel{k},fieldarray{k});
            end
        end
        s_out = currlevel{1};
        
    else
        s_out = rmfield(s,fieldarray(1));
    end

end

