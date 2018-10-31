function gait_gui_draw(hAx,hObject, eventdata, handles)

% Load the sysplotter configuration information
load sysplotter_config


% Bring selected axis to foreground (use gca if this is not specified)
if ~exist('hAx','var') || isempty(hAx)
    hAx = gca;
end

found_figure = 0;
hSearch = hAx;
while ~found_figure
    hSearch = get(hSearch,'Parent');
    if strcmpi(get(hSearch,'type'),'figure')
        found_figure = 1;
        hFig = hSearch;
    end
end
figure(hFig)
    
% Use the mouse to select a series of points
[alpha1,alpha2,button] = ginputc('AxHandle',hAx,'ShowPoints',true,'Color',Colorset.spot);

% Unless a non-primary button was used for the last click, copy the first
% point to the end to make a closed loop
if button(end) == 1
    alpha1 = [alpha1; alpha1(1)];
    alpha2 = [alpha2; alpha2(1)];
    
    % Set the period of the gait as 2pi
    period = 2*pi;
    
    % Evenly space the points along one period
    t = linspace(0,period,numel(alpha1));
    
    % Fit a periodic spline to the selected points; the endslopes are found
    % by averaging the positions of the points before and after the join
    endslope1 = (alpha1(2)-alpha1(end-1))/(t(end)-t(end-2));
    endslope2 = (alpha2(2)-alpha2(end-1))/(t(end)-t(end-2));
    spline_alpha1 = spline(t,[endslope1;alpha1(:);endslope1]);
    spline_alpha2 = spline(t,[endslope2;alpha2(:);endslope2]);
    
    
else
    
    % Set the period of the open motion as 1
    period = 1;
    
    % Evenly space the points along time from zero to 2pi
    t = linspace(0,period,numel(alpha1));
   
    % Fit a non-periodic spline to the selected points
    spline_alpha1 = spline(t,alpha1(:));
    spline_alpha2 = spline(t,alpha2(:));

end

% Upsample and plot to show gait to user
n_plot = 100;
t_plot = linspace(0,period,n_plot);

alpha1_plot = ppval(spline_alpha1,t_plot);
alpha2_plot = ppval(spline_alpha2,t_plot);

% y=optimalgaitgenerator(s,2,100);
% alpha1_optimalplot = y(:,1)';
% alpha2_optimalplot = y(:,2)';

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

gaitline = line('Parent',hAx,'XData',alpha1_plot,'YData',alpha2_plot,'ZData',maxZ*ones(size(alpha1_plot)),'Color',Colorset.spot,'LineWidth',5);

%%%% Ask the user for a filename
current_dir = pwd; % Remember where we started
cd(shchpath)       % Move to shape change directory

% Make sure that the user selects a valid filename and doesn't try putting
% the file somewhere other than the Shape Changes directory

while 1

    [paramfilename,pathname] = uiputfile('*.mat','Select a file name starting with ''params_''','params_.mat');
    
    if isequal(paramfilename,0) && isequal(pathname,0) %if both filename and pathname are 0, user hit cancel
        usercancel = 1;
        break
    else
        [~,paramfilenamebare,ext] = fileparts(paramfilename); % Break open the filename to get the extension
    end
    
    if ( strncmpi(paramfilename,'params_',7) && strcmpi(ext,'.mat') )...
                    && strcmpi(fullfile(pathname,paramfilename),fullfile(shchpath,paramfilename))
        usercancel = 0;
        break
    else
        disp('Must choose a filename starting with ''params_'' and in the Shape_Changes directory')
    end
        
end

cd(current_dir)    % Go back to original directory
%%%%

% If the user didn't hit cancel, save the data and create a shchf file that
% reads the data and interprets it as a gait.
if ~usercancel
    
    % Save the data to a parameters file
    save(fullfile(shchpath,paramfilename),'alpha1','alpha2','t')
    
    % Create the file if it doesn't already exist; future work could be
    % more fine-grained about this (making sure that any patches are
    % applied vs not overwriting any hand-edits the user made) and allowing
    % user to enter a prettier string for the display name here.
    
   
    shchfile = fullfile(shchpath,['shchf_' paramfilenamebare]);
    if ~exist([shchfile '.m'],'file') || ~isfunction(shchfile)
        
        gait_gui_draw_make_shchf(paramfilenamebare, paramfilenamebare,2)
        
    end
    
    refresh_handle = findall(0,'tag','refresh_gui');  % Get the handle for the button
    refresh_handle.Callback(refresh_handle,0)       % Push the refresh button
    
else
    
    return
    
end

shch_index = get(handles.shapechangemenu,'Value');
shch_names = get(handles.shapechangemenu,'UserData');

current_shch = shch_names{shch_index};

[rn,cn]=find(strcmp(shch_names,['shchf_' paramfilenamebare]));
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