function [present] = isfieldnested(s,fields)
%isfieldnested:
% Description: This function checks for the presenence of fields in a
% nested structure
% Inputs:
% s: structure being explored
% fields: cell array of fields to be checked

present = 1;

if isfield(s,fields{1})
    if(size(fields,2)>1)
        for k = 2:size(fields,2)
            %check if the fields are present
            if isfield(getfield(s,fields{1:k-1}),fields{k})
                
            else
                present = 0;
                break;
            end
        end
    end
else
    present = 0;
end


end

