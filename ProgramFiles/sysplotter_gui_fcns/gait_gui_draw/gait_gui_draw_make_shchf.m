function gait_gui_draw_make_shchf(paramfilename, displayname, n_dim, isfamily)

load sysplotter_config

if(nargin < 4)
    isfamily = 0;
end

% Open the template file
if (n_dim >= 2) && (n_dim <= 5)
    if (isfamily)
        gait_gui_draw_template_file = 'gait_family_gui_draw_template'...
            + string(n_dim) + '.txt';
    else
        gait_gui_draw_template_file = 'gait_gui_draw_template'...
            + string(n_dim) + '.txt';
    end
    fidi = fopen(fullfile(fileparts(which('gait_gui_draw')),...
        gait_gui_draw_template_file));
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