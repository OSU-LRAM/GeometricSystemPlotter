function v = verify_configdir(targetdir)

    % Verify that the target directory has the necessary subdirectories
    targetcontents = dir(targetdir);
    dirlist = {targetcontents.name};
    
    % If necessary subdirectories are in place, return true
    v = ~isempty(strcmp('Shape_Changes',dirlist)) && ...
            ~isempty(strcmp('Systems',dirlist));



end