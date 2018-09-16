function gait_gui_draw_make_shchf(paramfilename, displayname)

load sysplotter_config

% Open the template file
fidi = fopen(fullfile(sysplotterpath,'Utilities','gait_gui_draw',...
    'gait_gui_draw_template.txt'));

% Create the output file
fido = fopen(fullfile(shchpath,['shchf_' paramfilename '.m']),'w');


while ~feof(fidi)
    
	% Read the next line from the input file
	s = fgetl(fidi);

    % Insert the shchf filename
    s = strrep(s,'AA_SHCHFILENAME',['shchf_' paramfilename]);
    
    % Give this file the correct paramfile
    s = strrep(s,'AA_PARAMSNAME',paramfilename);
    
    % Insert the displayname
    s = strrep(s,'AA_DISPLAYNAME',displayname);
    
    % Copy put the line into the new file
    fprintf(fido,'%s\n',s);
    
end