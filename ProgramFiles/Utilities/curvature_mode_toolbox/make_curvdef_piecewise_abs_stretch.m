% test template creator for extensible serpenoid
% 5/16
% 7 PM - Huzzah! It works! 1000 sysf scripts generated and none of it by hand >:D 
% 5/25 edited fopen line to include sprintf so files are easy to tell apart
%      - ok I tried to edit it. Gives funny file names. Will make do with
%      original files for now. 
% 5/30 - Added lines to test just sections of b. should streamline
%      integration later. Make sure the correct lines are selected!
% 6/19 - Re-running for piecewise files. Make sure all the right files are
%      selected!
%      - Totally worked. Just had to not be an idiot with the wrong things
%      selected. Now I have 300 system files for the piecewise system :D

%% Make sure the correct options are selected! %%
function make_curvdef_piecewise_abs_stretch(b,sysplotterpath,syspath)

% b = 0.3;
% template file, in same directory
template = 'curvdef_stretch_piecewise_abs_template.txt'; % For continuous curvature piecewise
%           curvdef_stretch_piecewise_abs_template
% template = 'stretch_test_template.txt';      % For serpenoid 
% template = 'stretch_test_section_template.txt';

% curvdef_tempdir = fullfile(sysplotterpath,'Utilities','curvature_mode_toolbox',...
%     'curvdef_tempdir');
% 
% mkdir(curvdef_tempdir)

% generate a new file for each b value
for ind = 1:numel(b)
    
    % open both the new file and the template
    % fullfile(curvdef_tempdir,[fnamelist{idx} '.m'])
    fido = fopen(fullfile(syspath,['curv_piecewise_const_stretch_abs_b_',num2str(b(ind)*1000),'.m']),'w');
%                                       'piecewise_const_L_rms_b_',num2str(stretch_const*1000)
%     fido = fopen(fullfile(syspath,['curv_piecewise_const_L_rms_',num2str(b(ind)),'.m']),'w');
    fidi = fopen(fullfile(sysplotterpath,'Utilities','curvature_mode_toolbox',template));
    
    while ~feof(fidi) % while not at the end of the file read previously
        
        % read next line
        s = fgetl(fidi);
        
        % find locations of stretch_const in that line
        st_loc_d = strfind(s,'%d');
        st_loc_s = strfind(s,'%s');
        
        % replace "%d" or "%s" with values from b, either alone or *1000
        if ~isempty(st_loc_d)
                
            % find where the exact value of b is needed
            s = sprintf(s,b(ind));
            
        end
        
        if ~isempty(st_loc_s)
            
            % find where the *1000 value is needed 
%             s = sprintf(s,num2str(b(ind)*1000));
            s = sprintf(s,num2str(b(ind)*1000));
            
        end
        
        % write the next line in the file
        fprintf(fido,'%s\n',s);

        
    end
    
    % close the things
    fclose(fidi);
    fclose(fido);
end
end







