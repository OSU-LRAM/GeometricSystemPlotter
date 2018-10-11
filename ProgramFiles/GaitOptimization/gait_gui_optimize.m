function gait_gui_optimize(hAx,hObject,eventdata,handles)

% Load the sysplotter configuration information
load sysplotter_config

% Upsample and plot to show gait to user

system_index = get(handles.systemmenu,'Value');
system_names = get(handles.systemmenu,'UserData');

current_system = system_names{system_index};

shch_index = get(handles.shapechangemenu,'Value');
shch_names = get(handles.shapechangemenu,'UserData');

current_shch = shch_names{shch_index};

% insert warning/adjustment if file doesn't exist here.
% save(fullfile([sysplotter_inputpath '\Shape_Changes'],paramfilename),'alpha1','alpha2','t');

g1 = fullfile(shchpath,strcat(current_shch(7:end),'.mat'));
g2 = fullfile(datapath,strcat(current_system,'__',current_shch,'.mat'));
if exist(g1,'file') == 2
    load(g1);
elseif exist(g2,'file') == 2
    load(g2);
    if numel(p.phi_fun_full) == 1
        t = p.time_full{1,1};
        a12 = p.phi_fun_full{1,1}(t);
        alpha1 = a12(:,1);
        alpha2 = a12(:,2);
    else
        error('Selected gait file has more than one cycle. Optimizer only works on single cycles.')
    end
end

% Interpolating to increase number of points forming the gait
endslope1 = (alpha1(2)-alpha1(end-1))/(t(end)-t(end-2));
endslope2 = (alpha2(2)-alpha2(end-1))/(t(end)-t(end-2));
spline_alpha1 = spline(t,[endslope1;alpha1(:);endslope1]);
spline_alpha2 = spline(t,[endslope2;alpha2(:);endslope2]);
period = 2*pi;
n_plot = 100;
t_plot = linspace(0,period,n_plot+1);
alpha1_plot = ppval(spline_alpha1,t_plot);
alpha2_plot = ppval(spline_alpha2,t_plot);

% loading sysf file 
f=fullfile(datapath,strcat(current_system,'_calc.mat'));
load(f);

% Calling the optimizer
lb=0.95*[s.grid_range(1)*ones(n_plot+1,1);s.grid_range(3)*ones(n_plot+1,1)];%0.9 was points value
ub=0.95*[s.grid_range(2)*ones(n_plot+1,1);s.grid_range(4)*ones(n_plot+1,1)];
if strcmpi(s.system_type,'drag')
    % Use drag-dominated gait generator
    y=optimalgaitgenerator(s,2,n_plot,alpha1_plot,alpha2_plot,lb,ub,hAx,hObject,eventdata,handles);
elseif strcmpi(s.system_type,'inertial')
    % Use inertia-based gait generator
end
    
alpha1 = [y(1:100)',y(1)]';
alpha2 = [y(101:200)',y(101)]';
t=t_plot;
tnew=t(1):(t(2)-t(1))/4:t(end);
alpha12=interp1(t,alpha1,tnew,'spline');
alpha22=interp1(t,alpha2,tnew,'spline');


% Provide zdata to line if necessary
maxZ = 0;
hAxChildren = get(hAx,'Children');
if ~isempty(hAxChildren)
   for idx = 1:numel(hAxChildren)
       if ~isempty(hAxChildren(idx).ZData)
           maxZ = max(maxZ,max(hAxChildren(idx).ZData(:)));
       end
   end
end

gaitline = line('Parent',hAx,'XData',alpha12,'YData',alpha22,'ZData',maxZ*ones(size(alpha1)),'Color',Colorset.spot,'LineWidth',5);

%%%% Ask the user for a filename
current_dir = pwd; % Remember where we started
cd(shchpath)       % Move to shape change directory


% paramfilename=[current_shch(7:end) '.mat'];
% [~,paramfilenamebare,ext] = fileparts(paramfilename);

cd(current_dir)    % Go back to original directory
%%%%

sysf_func = str2func(current_system);
shch_func = str2func(current_shch);
paramfiledisplaytext = ['Opt: [',sysf_func('name'),'] [',shch_func('name'),'] Xeff'];
paramfiletext = ['opt_',current_system(6:end),'_',current_shch(7:end),'_Xeff'];

% If the user didn't hit cancel, save the data and create a shchf file that
% reads the data and interprets it as a gait.
%  if ~usercancel
    
    % Save the data to a parameters file
%     save(fullfile(shchpath,paramfilename),'alpha1','alpha2','t')
%     save(fullfile(shchpath,strcat(paramfilenamebare,'_optimal.mat')),'alpha1','alpha2','t')
    save(fullfile(shchpath,strcat(paramfiletext,'.mat')),'alpha1','alpha2','t')
    
    % Create the file if it doesn't already exist; future work could be
    % more fine-grained about this (making sure that any patches are
    % applied vs not overwriting any hand-edits the user made) and allowing
    % user to enter a prettier string for the display name here.
    
    % ['shchf_' paramfilenamebare '.m']
    if ~exist(fullfile(shchpath,['shchf_' paramfiletext '.m']),'file')
        % [paramfilenamebare '_optimal'], [paramfilenamebare '_optimal']
        gait_gui_draw_make_shchf(paramfiletext,paramfiledisplaytext)
        
    end
    
    refresh_handle = findall(0,'tag','refresh_gui');  % Get the handle for the button
    refresh_handle.Callback(refresh_handle,0)       % Push the refresh button
    
%  end
  
shch_index = get(handles.shapechangemenu,'Value');
shch_names = get(handles.shapechangemenu,'UserData');

current_shch = shch_names{shch_index};
%                              ['shchf_' paramfilenamebare '_optimal']
[rn,cn]=find(strcmp(shch_names,['shchf_' paramfiletext]));
set(handles.shapechangemenu,'Value',rn(1));
active=0;
% shapechangemenu_Callback(hObject, eventdata, handles,active)
enable_disable_shch_plots(hObject,eventdata,handles)

% Get the last plot pushbutton used
if isfield(handles.figure1.UserData,'lastpushbutton')
    lastpushbutton = handles.figure1.UserData.lastpushbutton;
else
    lastpushbutton = 'plotpushbutton1';
end

 plot_info = plotpushbutton_Callback(findall(0,'tag',lastpushbutton), eventdata, handles);   
    

end