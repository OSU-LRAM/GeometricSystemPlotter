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

g=fullfile(shchpath,strcat(current_shch(7:end),'.mat'));
load(g);

endslope1 = (alpha1(2)-alpha1(end-1))/(t(end)-t(end-2));
endslope2 = (alpha2(2)-alpha2(end-1))/(t(end)-t(end-2));
spline_alpha1 = spline(t,[endslope1;alpha1(:);endslope1]);
spline_alpha2 = spline(t,[endslope2;alpha2(:);endslope2]);
period = 2*pi;



n_plot = 100;
t_plot = linspace(0,period,n_plot+1);

alpha1_plot = ppval(spline_alpha1,t_plot);
alpha2_plot = ppval(spline_alpha2,t_plot);

f=fullfile(datapath,strcat(current_system,'_calc.mat'));
load(f);
lb=0.8*[s.grid_range(1)*ones(n_plot+1,1);s.grid_range(3)*ones(n_plot+1,1)];%0.9 was points value
ub=0.8*[s.grid_range(2)*ones(n_plot+1,1);s.grid_range(4)*ones(n_plot+1,1)];
y=optimalgaitgenerator(s,2,n_plot,alpha1_plot,alpha2_plot,lb,ub);
alpha1 = [y(1:100)',y(1)]';
alpha2 = [y(101:200)',y(101)]';
t=t_plot;

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

gaitline = line('Parent',hAx,'XData',alpha1,'YData',alpha2,'ZData',maxZ*ones(size(alpha1)),'Color',Colorset.spot,'LineWidth',5);

%%%% Ask the user for a filename
current_dir = pwd; % Remember where we started
cd(shchpath)       % Move to shape change directory


paramfilename=[current_shch(7:end) '.mat'];
[~,paramfilenamebare,ext] = fileparts(paramfilename);

cd(current_dir)    % Go back to original directory
%%%%

% If the user didn't hit cancel, save the data and create a shchf file that
% reads the data and interprets it as a gait.
%  if ~usercancel
    
    % Save the data to a parameters file
%     save(fullfile(shchpath,paramfilename),'alpha1','alpha2','t')
    save(fullfile(shchpath,strcat(paramfilenamebare,'_optimal.mat')),'alpha1','alpha2','t')
    
    % Create the file if it doesn't already exist; future work could be
    % more fine-grained about this (making sure that any patches are
    % applied vs not overwriting any hand-edits the user made) and allowing
    % user to enter a prettier string for the display name here.
    
    
    if ~exist(fullfile(shchpath,['shchf_' paramfilenamebare '_optimal.m']),'file')
        
        gait_gui_draw_make_shchf([paramfilenamebare '_optimal'], [paramfilenamebare '_optimal'] )
        
    end
    
    refresh_handle = findall(0,'tag','refresh_gui');  % Get the handle for the button
    refresh_handle.Callback(refresh_handle,0)       % Push the refresh button
    
%  end
  
 shch_index = get(handles.shapechangemenu,'Value');
shch_names = get(handles.shapechangemenu,'UserData');

current_shch = shch_names{shch_index};

[rn,cn]=find(strcmp(shch_names,['shchf_' paramfilenamebare '_optimal']));
set(handles.shapechangemenu,'Value',rn(1));
active=0;
% shapechangemenu_Callback(hObject, eventdata, handles,active)
enable_disable_shch_plots(hObject,eventdata,handles)
 plot_info = plotpushbutton_Callback(findall(0,'tag','plotpushbutton3'), eventdata, handles);   
    

end