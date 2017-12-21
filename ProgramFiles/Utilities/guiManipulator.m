% function guiManipulator(b_array)

% test gui manipulation

% 5/26 - tentative success :D going to run with 20 trials beginning at
%   b=0.1. Probably will take a while. Future note; should make this a function, 
%   so I can input what kind of data I want calculated. 
%   note: ran 31 systems in 526 seconds, with one that didn't update. Looks
%   like it takes ~17 seconds to update one system file. 
%   Animation test confirms the data is what I expected.

% 6/19 - Added for loop in data thing to allow searching for a specific
%   system within sysplotter's system list. Sysplotter still needs the
%   right shape change selected, but at least this makes the system
%   selection a little easier. 
%   - Attempted to make this into a function. Apparently Matlab doesn't
%   like loading variables into a static workspace. Must find a way around
%   that. Current workaround is to have a pre-existing b_array already
%   defined in the workspace. 

% clearvars -except b_array sysplotter_inputpath

if ~exist('b_array','var')
   b_dummy = 0:0.01:0.1;
%     b_dummy = round(linspace(0,0.2,10),3);
else
   b_dummy = b_array; 
end
userhandles = guidata(sysplotter); % sysplotter should be open
% load('C:\Users\Jacquelin\OneDrive - Oregon State University\OSU\0 Research Other\Matlab\0 GSP files\guiManipulation\3\hObject.mat');
% 
% %extract name of source
% source_name = get(hObject,'Tag');

%extract the column number as text
% source_number_text = source_name(end);

%get the checkbox values
% [box_names, box_active, box_values, box_enabled, plot_types, ...
%     plot_subtypes,merged_plot_subtypes, plot_style] =...
%     get_box_values(source_number_text,userhandles);

%Get the system and shape change info
system_index = get(userhandles.systemmenu,'Value');
system_names = get(userhandles.systemmenu,'UserData');

current_system = system_names{system_index};
% current_system = strfind(userhandles.systemmenu.String,'Serpenoid Stretchable Snake, b = 0');

shch_index = get(userhandles.shapechangemenu,'Value');
shch_names = get(userhandles.shapechangemenu,'UserData');

current_shch = shch_names{shch_index};

% Get the desired vector and scalar plotting resolutions
resolution.vector = str2num(get(userhandles.vectorresolution,'String'));
resolution.scalar = str2num(get(userhandles.scalarresolution,'String'));
resolution.vector_range = get(userhandles.vectorresolution,'UserData');
resolution.scalar_range = get(userhandles.scalarresolution,'UserData');

disp(['Shape change: ',current_shch])
disp('Updating system data:')

tic
for i = 1:numel(b_dummy)
    if strcmp(length_sys,'original')
        if strcmp(snake_sys,'serpenoid')
            if strcmp(lambda_func,'quad')
                current_system = ['sysf_serpenoid_extendable_',num2str(b_dummy(i)*1000)];
            elseif strcmp(lambda_func,'abs')
                current_system = ['sysf_serpenoid_extendable_abs_',num2str(b_dummy(i)*1000)];
            end
        elseif strcmp(snake_sys,'piecewise')
            if strcmp(lambda_func,'quad')
                current_system = ['sysf_piecewise_const_stretch_b_',num2str(b_dummy(i)*1000)];
            elseif strcmp(lambda_func,'abs')
                current_system = ['sysf_piecewise_const_stretch_abs_b_',num2str(b_dummy(i)*1000)];
            end
        end
    elseif strcmp(length_sys,'Lrms')
        if strcmp(snake_sys,'serpenoid')
            if strcmp(lambda_func,'quad')
                current_system = ['sysf_serpenoid_quad_Lrms_',num2str(b_dummy(i)*1000)];
            elseif strcmp(lambda_func,'abs')
                current_system = ['sysf_serpenoid_abs_Lrms_',num2str(b_dummy(i)*1000)];
            end
        elseif strcmp(snake_sys,'piecewise')
            if strcmp(lambda_func,'quad')
                current_system = ['sysf_piecewise_quad_Lrms_',num2str(b_dummy(i)*1000)];
            elseif strcmp(lambda_func,'abs')
                current_system = ['sysf_piecewise_abs_Lrms_',num2str(b_dummy(i)*1000)];
            end
        end
    end
    
    update_no_plot(current_system,current_shch);

%     clear ind
%     for j = 1:numel(system_names)
%         if isempty(strfind(system_names{j},current_system))
%             ind(j) = 0;
%         else
%             ind(j) = 1;
%         end
%     end
%     system_index = find(ind);
    
    %     system_index = system_index+1;
%     current_system = system_names{system_index};
    
end
t_end = toc; 

disp(['Finished refreshing data on ',num2str(i),' systems in ',num2str(t_end),' seconds'])

function update_no_plot(sys,shch) % add userhandles back in?

% declare the data directory
configfile = './sysplotter_config';
load(configfile);

% something here isn't letting the plots updata data properly.
%
% 5/26 dummy me forgot to copy the code that actually updates the data :P 

disp(['    ',sys,'...']);

%convert empty matrix input to 'null' string
if isempty(shch)
    
    shch = 'null';
    
end

% Determine which components need to be updated
if exist('handles','var')
    [update] = refresh_runinfo_Callback([], [], handles);
else
    [update] = decide_components(sys,shch);
end

% Recalculate the details if necessary

% Load a structure with the path names
pathnames = load(configfile);

if update.sys_init
    disp('        Initializing system');
    s = absolute_feval(fullfile(syspath, sys),'initialize',pathnames); %#ok<NASGU>
    save(fullfile(datapath, sys),'s')
end

if update.shch_init
    disp('        Initializing shape change');    
    p = absolute_feval(fullfile(shchpath, shch),'initialize',pathnames); %#ok<NASGU>
    save(fullfile(datapath, shch),'p')
end

if update.sys_calc
    disp('        Updating calculations');
    sys_calcsystem('calculate',sys);
end

if update.shch_calc
    disp('        Updating shape change');
    sys_calcpath('calculate',sys,shch);
end

disp('    done.');

end



% end


