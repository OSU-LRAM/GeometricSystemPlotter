% optimizerManipulator

% 9/28 - Utilizing Suresh's gait optimizer. I would like to have something
%   similar to the original manipulator, where I call the function Callback
%   and it runs on all systems. In this case, it would need to use the data
%   from the last to initialize the next. Or I could have it start with the
%   circle and continue to run, but that may take quite a long time. Who
%   knows, I'll just jump in.
% 10/7 - Started backbone for code. Discovered the findobj is pretty much
%   unnecessary, just use guidata(sysplotter) instead. Still need to add in
%   creating dummy gaits and actually calling the optimizer. (also I
%   totally discovered that a long time ago and just plain forgot :P )
% 10/9 - Mostly finished. Added in the optimizer callback, and added
%   options for which backbone and lambda function to enable so I don't
%   have to keep commenting and uncommenting like a scrub. Planning on
%   running it four times on 10 systems of each combination.
% 10/27 - Modified code to copy old optimal path into the new seed file, so
%   ideally the optimizer won't get caught up and give odd values.
% 11/13 - Note for future: turns out fminconn is doing something odd with
%   step sizes, so Suresh will have to write his own optimizer. Might have
%   to make changes to my code when that's complete.


clearvars -except sysplotter_inputpath

snake_sys = 'serpenoid';    % Which backbone to use?
% snake_sys = 'piecewise';

lambda_func = 'quad';       % Which lambda input function?
% lambda_func = 'abs';

% length_sys = 'Lrms';        % Lrms length or changeable?
length_sys = 'original';

save_figs = 0;              % should sysplotter save figures?

run_optimize = 1;           % should we run the optimizer on all systems?

handles = guidata(sysplotter);

seed_method = 'prev_opt';
% seed_method = 'circle'; 

% sysplotter('plotpushbutton_Callback',handles.plotpushbutton3,[],handles) % not this, but same format.

% Generate temp gait and data as an optimizer seed for first run
th = linspace(2*pi,0,101);
amp = 4.5;
alpha1 = amp*cos(th);
alpha2 = amp*sin(th);
t = flipud(th');

% b = 0:0.01:0.1;
if strcmp(lambda_func,'quad')
    b = round(linspace(0,0.101,21),3); % quadratic value b range
%     b = b(9:11);
elseif strcmp(lambda_func,'abs')
    b = round(linspace(0,0.318,21),3); % absolute value b range
    %         b = b(1);
end

% make sure the correct checkbox is checked
set(handles.CCFXoptcheckbox3,'Value',1)

% update null shch values
set(handles.systemmenu,'Value',1) % set no system for guiManipulator
set(handles.shapechangemenu,'Value',1) % set no shapechange

b_array = b;
guiManipulator;

% pause(1)
if run_optimize
    beep %#ok<*UNRCH>
    
%     old_param_file = 'piecewise_seed_b_40_optimal';
    
    for i = 1:numel(b)
        
        
        % savename = ['optimizer_circle_start_',num2str(b*1000),'.mat'];
        
        % figure out which file and system we're using
        
        if strcmp(length_sys,'original')
            if strcmp(snake_sys,'serpenoid') % serpenoid system
                if strcmp(lambda_func,'quad') % quadratic lambda function
                    paramfilename = ['serpenoid_seed_b_',num2str(b(i)*1000)];
                    current_system = ['sysf_serpenoid_extendable_',num2str(b(i)*1000)];
                elseif strcmp(lambda_func,'abs') % absolute lambda function
                    paramfilename = ['serpenoid_seed_abs_b_',num2str(b(i)*1000)];
                    current_system = ['sysf_serpenoid_extendable_abs_',num2str(b(i)*1000)];
                end
            elseif strcmp(snake_sys,'piecewise') % piecewise system
                if strcmp(lambda_func,'quad')
                    paramfilename = ['piecewise_seed_b_',num2str(b(i)*1000)];
                    current_system = ['sysf_piecewise_const_stretch_b_',num2str(b(i)*1000)];
                elseif strcmp(lambda_func,'abs')
                    paramfilename = ['piecewise_seed_abs_b_',num2str(b(i)*1000)];
                    current_system = ['sysf_piecewise_const_stretch_abs_b_',num2str(b(i)*1000)];
                end
            end
        elseif strcmp(length_sys,'Lrms')
            if strcmp(snake_sys,'serpenoid') % serpenoid system
                if strcmp(lambda_func,'quad') % quadratic lambda function
                    paramfilename = ['serpenoid_seed_Lrms_b_',num2str(b(i)*1000)];
                    current_system = ['sysf_serpenoid_quad_Lrms_',num2str(b(i)*1000)];
                elseif strcmp(lambda_func,'abs') % absolute lambda function
                    paramfilename = ['serpenoid_seed_abs_Lrms_b_',num2str(b(i)*1000)];
                    current_system = ['sysf_serpenoid_abs_Lrms_',num2str(b(i)*1000)];
                end
            elseif strcmp(snake_sys,'piecewise') % piecewise system
                if strcmp(lambda_func,'quad')
                    paramfilename = ['piecewise_seed_Lrms_b_',num2str(b(i)*1000)];
                    current_system = ['sysf_piecewise_quad_Lrms_',num2str(b(i)*1000)];
                elseif strcmp(lambda_func,'abs')
                    paramfilename = ['piecewise_seed_abs_Lrms_b_',num2str(b(i)*1000)];
                    current_system = ['sysf_piecewise_abs_Lrms_',num2str(b(i)*1000)];
                end
            end
        end
        
        disp(['system: ',current_system])
        disp(['parameter file: ',paramfilename])
        
        
        if strcmp(seed_method,'prev_opt') % save .mat file and create new shch file
            if i == 1
                save(fullfile([sysplotter_inputpath '\Shape_Changes'],paramfilename),'alpha1','alpha2','t');
            else
                old_optimal = [fullfile([sysplotter_inputpath '\Shape_Changes'],old_param_file),'.mat'];
                new_seed = [fullfile([sysplotter_inputpath '\Shape_Changes'],paramfilename),'.mat'];
                copyfile(old_optimal,new_seed)
            end
        elseif strcmp(seed_method,'circle') % use circles as seeds
            save(fullfile([sysplotter_inputpath '\Shape_Changes'],paramfilename),'alpha1','alpha2','t');
        end
        
        
        if strcmp(lambda_func,'quad')
            displayname = [snake_sys,' optimizer start circle, b = ',num2str(b(i))];
        elseif strcmp(lambda_func,'abs')
            displayname = [snake_sys,' optimizer start circle, b = ',num2str(b(i)),', abs lambda func'];
        end
        
        
        disp('    creating shape change file')
        gait_gui_draw_make_shchf(paramfilename, displayname);
        
        % refresh sysplotter so it can see the new shch
        disp('    refreshing sysplotter')
        sysplotter('refresh_gui_Callback',handles.refresh_gui,[],handles)
        
        %     disp('    plotting null shch')
        %     set(handles.shapechangemenu,'Value',1)
        
        system_names = get(handles.systemmenu,'UserData');
        shch_names = get(handles.shapechangemenu,'UserData');
        
        disp('    setting current system and shch')
        % set current system
        dummy_idx_start = strfind(system_names,'end_recent');
        idx_start = find(not(cellfun('isempty',dummy_idx_start)));
        dummy_idx = strfind(system_names(idx_start+1:end),current_system);
        %     dummy_idx = dummy_idx(idx_start+1:end);
        idx = find(not(cellfun('isempty',dummy_idx)));
        if isempty(idx)
            error(['System ',current_system,' not found.'])
        end
        set(handles.systemmenu,'Value',idx(1)+idx_start)
        % set current shapechange to optimizer seed
        dummy_idy_start = strfind(shch_names,'end_recent');
        idy_start = find(not(cellfun('isempty',dummy_idy_start)));
        dummy_idy = strfind(shch_names(idy_start+1:end),['shchf_',paramfilename]);
        %     dummy_idy = dummy_idy(idy_start+1:end);
        idy = find(not(cellfun('isempty',dummy_idy)));
        if isempty(idy)
            error(['parameter File ',paramfilename,' not found.'])
        end
        set(handles.shapechangemenu,'Value',idy(1)+idy_start)
        
        % refresh again so the plot buttons are enabled?
        % disp('    refreshing sysplotter')
        % sysplotter('refresh_gui_Callback',handles.refresh_gui,[],handles)
        
        % set plot button to enabled just in case it isn't
        disp('    plotting new shch')
        set(handles.plotpushbutton3,'enable','on')
        
        % because it errors sometimes
        set_plot_resolution_from_file(handles)
        
        % Plot the new shch so the optimizer can see it and run
        set(handles.CCFXoptcheckbox3,'Value',1)
        sysplotter('plotpushbutton_Callback',handles.plotpushbutton3,[],handles)
        % set(handles.CCFXoptcheckbox3,'Value',0)
        
        % call the optimizer
        disp('    running optimizer: prepare for huge block of text')
        tic
        sysplotter('OptimizeButton_Callback',handles.OptimizeButton,[],handles)
        opt_time(i) = toc; %#ok<SAGROW>
        
        % clear original seed to save filespace
        % delete(fullfile([sysplotter_inputpath '\Shape_Changes'],[paramfilename,'.mat']));
        % delete(fullfile([sysplotter_inputpath '\Shape_Changes'],['shchf_',paramfilename,'.m']));
        
        % disp(['finished optimizing path for ',snake_sys,' backbone with b = ',num2str(b),' in ',num2str(toc),' seconds.'])
        
        % save the old parameter file as the new seed for iterations after i = 1
        old_param_file = [paramfilename,'_optimal'];
        
    end
    
    save(fullfile([sysplotter_inputpath '\sysplotter_data'],[snake_sys,'_',lambda_func,'_',length_sys,'_runtimes.mat']),'opt_time');
    disp(['Finished running ',snake_sys,' ',lambda_func,' ',length_sys,' with seed method ',seed_method])
    clearvars -except sysplotter_inputpath opt_time
else
    disp('Not optimizing.')
end

beep