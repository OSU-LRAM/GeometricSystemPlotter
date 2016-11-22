%%%%%%%%%%%%%%
%Get date function
function sdate = getfiledate(dirstruct, name)

    %index of file getting date for
    fileindex = strmatch(name, {dirstruct.name},'exact');
    
    %if the file exists, get its date
    if ~isempty(fileindex)
        
        sdate = dirstruct(fileindex).datenum;
        
    else
        
        sdate = 0;
        
    end
    
end
