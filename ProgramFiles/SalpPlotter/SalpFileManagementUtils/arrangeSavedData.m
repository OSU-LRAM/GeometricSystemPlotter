%Rearranges a cell array based on last input for 'recently used' dropdowns
function [nameList,fileList] = arrangeSavedData(nameSaved,fileSaved,oldNameList,oldFileList,nSave)

    %If cell array is too long delete elements until it's not
    while numel(oldNameList) > nSave    
        oldNameList(nSave) = [];
        oldFileList(nSave) = [];
    end

    %If we're saving something that's already in the recent list
    if any(strcmpi(nameSaved,oldNameList))

        nameList = oldNameList;
        fileList = oldFileList;

        %Delete where it currently appears
        appearIndex = find(strcmpi(nameSaved,oldNameList),1);
        nameList(appearIndex) = [];
        fileList(appearIndex) = [];

        %And put it at the top of the list
        nameList = concatenateCellArrays({nameSaved},nameList);
        fileList = concatenateCellArrays({fileSaved},fileList);

    %If what we're saving is not currently in the list
    else
        
        %If the list is going to be too long delete the last element
        if numel(oldNameList) == nSave    
            oldNameList(nSave) = [];
            oldFileList(nSave) = [];
        end

        %Add new value to the top of the list
        nameList = concatenateCellArrays({nameSaved},oldNameList);
        fileList = concatenateCellArrays({fileSaved},oldFileList);

    end

end
            
%Concatenates two cell vectors of strings
function combinedCellArray = concatenateCellArrays(c1,c2)
    combinedCellArray = {};
    for i = 1:numel(c1)
        combinedCellArray{end+1} = c1{i};
    end
    for i = 1:numel(c2)
        combinedCellArray{end+1} = c2{i};
    end
end
