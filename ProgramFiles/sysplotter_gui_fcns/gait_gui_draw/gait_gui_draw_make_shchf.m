function gait_gui_draw_make_shchf(paramfilename, displayname,n_dim)

load sysplotter_config

% Open the template file
if n_dim==2
    fidi = fopen(fullfile(fileparts(which('gait_gui_draw')),...
    'gait_gui_draw_template2.txt'));
elseif n_dim==3
     fidi = fopen(fullfile(fileparts(which('gait_gui_draw')),...
    'gait_gui_draw_template3.txt'));
elseif n_dim==4
     fidi = fopen(fullfile(fileparts(which('gait_gui_draw')),...
    'gait_gui_draw_template4.txt'));
elseif n_dim==5
    fidi = fopen(fullfile(fileparts(which('gait_gui_draw')),...
    'gait_gui_draw_template5.txt'));
else
    error('Trying to make an shchf with an unsupported number of dimensions')
end

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