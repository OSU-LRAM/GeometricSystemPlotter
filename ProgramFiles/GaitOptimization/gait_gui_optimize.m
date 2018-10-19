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
if exist(g2,'file') == 2
    load(g2);
    n_dim=length(s.vecfield.eval.content.Avec_optimized(1,:));
    if numel(p.phi_fun_full) == 1
        t = p.time_full{1,1};
        a12 = p.phi_locus_full{1,1}.shape;
        if n_dim>=2
            alpha1 = a12(:,1);
            alpha2 = a12(:,2);
        end
        if n_dim>=3
            alpha3 = a12(:,3);
        end
        if n_dim>=4
            alpha4 = a12(:,4);
        end       
elseif exist(g1,'file') == 2
    load(g1);
    if exist('alpha4')==1
        n_dim=4;
    elseif exist('alpha3')==1
        n_dim=3;
    else
        n_dim=2;
    end
    else
        error('Selected gait file has more than one cycle. Optimizer only works on single cycles.')
    end
end

if n_dim==2
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
elseif n_dim==3
    endslope1 = (alpha1(2)-alpha1(end-1))/(t(end)-t(end-2));
    endslope2 = (alpha2(2)-alpha2(end-1))/(t(end)-t(end-2));
    endslope3 = (alpha3(2)-alpha3(end-1))/(t(end)-t(end-2));
    spline_alpha1 = spline(t,[endslope1;alpha1(:);endslope1]);
    spline_alpha2 = spline(t,[endslope2;alpha2(:);endslope2]);
    spline_alpha3 = spline(t,[endslope3;alpha3(:);endslope3]);
    period = 2*pi;
    n_plot = 100;
    t_plot = linspace(0,period,n_plot+1);
    alpha1_plot = ppval(spline_alpha1,t_plot);
    alpha2_plot = ppval(spline_alpha2,t_plot);
    alpha3_plot = ppval(spline_alpha3,t_plot);
elseif n_dim==4
    endslope1 = (alpha1(2)-alpha1(end-1))/(t(end)-t(end-2));
    endslope2 = (alpha2(2)-alpha2(end-1))/(t(end)-t(end-2));
    endslope3 = (alpha3(2)-alpha3(end-1))/(t(end)-t(end-2));
    endslope4 = (alpha4(2)-alpha4(end-1))/(t(end)-t(end-2));
    spline_alpha1 = spline(t,[endslope1;alpha1(:);endslope1]);
    spline_alpha2 = spline(t,[endslope2;alpha2(:);endslope2]);
    spline_alpha3 = spline(t,[endslope3;alpha3(:);endslope3]);
    spline_alpha4 = spline(t,[endslope4;alpha4(:);endslope3]);
    period = 2*pi;
    n_plot = 100;
    t_plot = linspace(0,period,n_plot+1);
    alpha1_plot = ppval(spline_alpha1,t_plot);
    alpha2_plot = ppval(spline_alpha2,t_plot);
    alpha3_plot = ppval(spline_alpha3,t_plot);
    alpha4_plot = ppval(spline_alpha4,t_plot);
end
% loading sysf file 
f=fullfile(datapath,strcat(current_system,'_calc.mat'));
load(f);

if n_dim==2
    % Calling the optimizer
    lb=0.95*[s.grid_range(1)*ones(n_plot+1,1);s.grid_range(3)*ones(n_plot+1,1)];%0.9 was points value
    ub=0.95*[s.grid_range(2)*ones(n_plot+1,1);s.grid_range(4)*ones(n_plot+1,1)];
    y=optimalgaitgenerator(s,2,n_plot,alpha1_plot,alpha2_plot,lb,ub);
    alpha1 = [y(1:100)',y(1)]';
    alpha2 = [y(101:200)',y(101)]';
    t=t_plot;
    tnew=t(1):(t(2)-t(1))/4:t(end);
    alpha12=interp1(t,alpha1,tnew,'spline');
    alpha22=interp1(t,alpha2,tnew,'spline');
elseif n_dim==3
    % Calling the optimizer
    lb=0.95*[s.grid_range(1)*ones(n_plot+1,1);s.grid_range(3)*ones(n_plot+1,1);s.grid_range(5)*ones(n_plot+1,1)];%0.9 was points value
    ub=0.95*[s.grid_range(2)*ones(n_plot+1,1);s.grid_range(4)*ones(n_plot+1,1);s.grid_range(6)*ones(n_plot+1,1)];
    y=optimalgaitgenerator3(s,3,n_plot,alpha1_plot,alpha2_plot,alpha3_plot,lb,ub);
    alpha1 = [y(1:100)',y(1)]';
    alpha2 = [y(101:200)',y(101)]';
    alpha3 = [y(201:300)',y(201)]';
    t=t_plot;
    tnew=t(1):(t(2)-t(1))/4:t(end);
    alpha12=interp1(t,alpha1,tnew,'spline');
    alpha22=interp1(t,alpha2,tnew,'spline');
elseif n_dim==4
    % Calling the optimizer
    lb=0.95*[s.grid_range(1)*ones(n_plot+1,1);s.grid_range(3)*ones(n_plot+1,1);s.grid_range(5)*ones(n_plot+1,1);s.grid_range(7)*ones(n_plot+1,1)];%0.9 was points value
    ub=0.95*[s.grid_range(2)*ones(n_plot+1,1);s.grid_range(4)*ones(n_plot+1,1);s.grid_range(6)*ones(n_plot+1,1);s.grid_range(8)*ones(n_plot+1,1)];
    y=optimalgaitgenerator4(s,4,n_plot,alpha1_plot,alpha2_plot,alpha3_plot,alpha4_plot,lb,ub);
    alpha1 = [y(1:100)',y(1)]';
    alpha2 = [y(101:200)',y(101)]';
    alpha3 = [y(201:300)',y(201)]';
    alpha4 = [y(301:400)',y(301)]';
    t=t_plot;
    tnew=t(1):(t(2)-t(1))/4:t(end);
    alpha12=interp1(t,alpha1,tnew,'spline');
    alpha22=interp1(t,alpha2,tnew,'spline');
end


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

gaitline = line('Parent',hAx,'XData',alpha12,'YData',alpha22,'ZData',maxZ*ones(size(alpha12)),'Color',Colorset.spot,'LineWidth',5);

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
if n_dim==2
    save(fullfile(shchpath,strcat(paramfiletext,'.mat')),'alpha1','alpha2','t')
elseif n_dim==3
    save(fullfile(shchpath,strcat(paramfiletext,'.mat')),'alpha1','alpha2','alpha3','t')
elseif n_dim==4
    save(fullfile(shchpath,strcat(paramfiletext,'.mat')),'alpha1','alpha2','alpha3','alpha4','t')
end
    
    % Create the file if it doesn't already exist; future work could be
    % more fine-grained about this (making sure that any patches are
    % applied vs not overwriting any hand-edits the user made) and allowing
    % user to enter a prettier string for the display name here.
    
    % ['shchf_' paramfilenamebare '.m']
    if ~exist(fullfile(shchpath,['shchf_' paramfiletext '.m']),'file')
        % [paramfilenamebare '_optimal'], [paramfilenamebare '_optimal']
        gait_gui_draw_make_shchf(paramfiletext,paramfiledisplaytext,n_dim)
        
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