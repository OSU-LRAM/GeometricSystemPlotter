function gait_gui_optimize(hAx,hObject,eventdata,handles,OptFamily)

% Load the sysplotter configuration information
load sysplotter_config

% if the function does not take in isStepOptimzer, set it as 0.
if(nargin < 5)
    OptFamily = 0;
end

optimization_handles = handles;
handles = handles.main_gui;

% Upsample and plot to show gait to user

system_index = get(handles.systemmenu,'Value');
system_names = get(handles.systemmenu,'UserData');

current_system = system_names{system_index};

shch_index = get(handles.shapechangemenu,'Value');
shch_names = get(handles.shapechangemenu,'UserData');

current_shch = shch_names{shch_index};

stretch = get(handles.stretchmenu,'Value')-1;

% insert warning/adjustment if file doesn't exist here.
% save(fullfile([sysplotter_inputpath '\Shape_Changes'],paramfilename),'alpha1','alpha2','t');

% Path to the shape change definition
g1 = fullfile(shchpath,strcat(current_shch(7:end),'.mat'));

% Path to the output of running the shape change through the system
g2 = fullfile(datapath,strcat(current_system,'__',current_shch,'.mat'));

% If the output of running the gait exists
if exist(g2,'file') == 2
    
    % Load the path data from the output file
    load(g2,'p','s');
    
    % Get the number of dimensions from the size of the system's shape
    % space
    n_dim=length(s.vecfield.eval.content.Avec_optimized(1,:));
    
    
    if numel(p.phi_fun_full) == 1
        t = p.time_full{1,1};
        alpha = p.phi_locus_full{1,1}.shape;
%             if n_dim>=2
%                 alpha1 = a12(:,1);
%                 alpha2 = a12(:,2);
%             end
%             if n_dim>=3
%                 alpha3 = a12(:,3);
%             end
%             if n_dim>=4
%                 alpha4 = a12(:,4);
%             end
    else
        error('Selected gait file has more than one cycle. Optimizer only works on single cycles.')
    end
else
    error('Execute the gait on the system before optimizing')
end

% % If the output of running the gait does not exist, but the 
% elseif exist(g1,'file') == 2
%     load(g1);
%     if exist('alpha4')==1
%         n_dim=4;
%     elseif exist('alpha3')==1
%         n_dim=3;
%     else
%         n_dim=2;
%     end
% else
%     error('Selected gait file has more than one cycle. Optimizer only works on single cycles.')
% end

% Interpolate gait to increase number of points forming the gait
endslope = (alpha(2,:)-alpha(end-1,:))/(t(end)-t(end-2));
spline_alpha = spline(t,[endslope;alpha;endslope]');
period = 2*pi;
n_plot = 100;
t_plot = linspace(0,period,n_plot+1);
alpha_plot = ppval(spline_alpha,t_plot)';



% if n_dim==2
%     % Interpolating to increase number of points forming the gait
%     endslope1 = (alpha1(2)-alpha1(end-1))/(t(end)-t(end-2));
%     endslope2 = (alpha2(2)-alpha2(end-1))/(t(end)-t(end-2));
%     spline_alpha1 = spline(t,[endslope1;alpha1(:);endslope1]);
%     spline_alpha2 = spline(t,[endslope2;alpha2(:);endslope2]);
%     period = 2*pi;
%     n_plot = 100;
%     t_plot = linspace(0,period,n_plot+1);
%     alpha1_plot = ppval(spline_alpha1,t_plot);
%     alpha2_plot = ppval(spline_alpha2,t_plot);
% elseif n_dim==3
%     endslope1 = (alpha1(2)-alpha1(end-1))/(t(end)-t(end-2));
%     endslope2 = (alpha2(2)-alpha2(end-1))/(t(end)-t(end-2));
%     endslope3 = (alpha3(2)-alpha3(end-1))/(t(end)-t(end-2));
%     spline_alpha1 = spline(t,[endslope1;alpha1(:);endslope1]);
%     spline_alpha2 = spline(t,[endslope2;alpha2(:);endslope2]);
%     spline_alpha3 = spline(t,[endslope3;alpha3(:);endslope3]);
%     period = 2*pi;
%     n_plot = 100;
%     t_plot = linspace(0,period,n_plot+1);
%     alpha1_plot = ppval(spline_alpha1,t_plot);
%     alpha2_plot = ppval(spline_alpha2,t_plot);
%     alpha3_plot = ppval(spline_alpha3,t_plot);
% elseif n_dim==4
%     endslope1 = (alpha1(2)-alpha1(end-1))/(t(end)-t(end-2));
%     endslope2 = (alpha2(2)-alpha2(end-1))/(t(end)-t(end-2));
%     endslope3 = (alpha3(2)-alpha3(end-1))/(t(end)-t(end-2));
%     endslope4 = (alpha4(2)-alpha4(end-1))/(t(end)-t(end-2));
%     spline_alpha1 = spline(t,[endslope1;alpha1(:);endslope1]);
%     spline_alpha2 = spline(t,[endslope2;alpha2(:);endslope2]);
%     spline_alpha3 = spline(t,[endslope3;alpha3(:);endslope3]);
%     spline_alpha4 = spline(t,[endslope4;alpha4(:);endslope3]);
%     period = 2*pi;
%     n_plot = 100;
%     t_plot = linspace(0,period,n_plot+1);
%     alpha1_plot = ppval(spline_alpha1,t_plot);
%     alpha2_plot = ppval(spline_alpha2,t_plot);
%     alpha3_plot = ppval(spline_alpha3,t_plot);
%     alpha4_plot = ppval(spline_alpha4,t_plot);
% end

% Direction to optimize: 1-x, 2-y, 3-theta
% Both buttons in Optimizer and Step-Optimizer panel have same tags
% They are distinguished by the array determined by the order of panel.
direction = 0;
if optimization_handles.xbutton(OptFamily+1).Value
    direction = 1;
elseif optimization_handles.ybutton(OptFamily+1).Value
    direction = 2;
elseif optimization_handles.thetabutton(OptFamily+1).Value
    direction = 3;
elseif (optimization_handles.steeringbutton.Value) && (OptFamily)
    direction = 4;
end

%Cost function to optimize over
costfunction = 'none';
if optimization_handles.pathlengthmetricbutton.Value
    costfunction = 'pathlength metric';
elseif optimization_handles.torquebutton.Value
    costfunction = 'torque';
elseif optimization_handles.covaccbutton.Value
    costfunction = 'covariant acceleration';
elseif optimization_handles.pathlengthcoordinatebutton.Value
    costfunction = 'pathlength coord';
elseif optimization_handles.accelerationbutton.Value
    costfunction = 'acceleration coord';
elseif optimization_handles.pathlengthmetric2button.Value
    costfunction = 'pathlength metric2';
elseif optimization_handles.powerqualitybutton.Value
    costfunction = 'power quality';
end

% Active constraint to optimzie
constraint = zeros(2,1);
if optimization_handles.otherdirectioncheckbox.Value
    constraint(1) = 1;
end
if optimization_handles.rotationcheckbox.Value
    constraint(2) = 1;
end

%Sanity check radiobutton selection
if strcmpi(costfunction,'none')
    error('Something went wrong with the costfunction selection');
end

% loading sysf file 
f=fullfile(datapath,strcat(current_system,'_calc.mat'));
load(f,'s');

%%%% Set lower and upper bounds on optimization

% Get the lower and upper bounds on the grid
lvals = s.grid_range(1:2:end);
uvals = s.grid_range(2:2:end);


% Get the center of the grid
mvals = (uvals + lvals)/2;

% Get the lower and upper grid bounds relative to the center of the grid
dlvals = lvals-mvals;
duvals = uvals-mvals;

% scale the differences by .95
sdlvals = dlvals * 0.95;
sduvals = duvals * 0.95;

% Add the scaled differences to the center
slvals = mvals+sdlvals;
suvals = mvals+sduvals;

% Tile these out to a matrix with as many rows as there are time points
lb = 0.95 * repmat(slvals,[(n_plot+1),1]);
ub = 0.95 * repmat(suvals,[(n_plot+1),1]);

% Stack the lower bound values and the upper bound values
lb = lb(:);
ub = ub(:);

%%%%% Call the optimizer
y = optimalgaitgenerator(s,n_dim,n_plot,alpha_plot,lb,ub,stretch,direction,costfunction,constraint,OptFamily,handles);

t = t_plot;
t_new = t(1):(t(2)-t(1))/4:t(end);

if(OptFamily)

    alpha_out = cell(4,1);
    for i = 1:4
        % reshape the output and add the start point to the end to close the loop
        alpha_out{i} = reshape(y{i},[numel(y{i})/n_dim,n_dim]);
        alpha_out{i} = [alpha_out{i};alpha_out{i}(1,:)];
    end
else
    % reshape the output and add the start point to the end to close the loop
    alpha_out = reshape(y,[numel(y)/n_dim,n_dim]);
    alpha_out = [alpha_out;alpha_out(1,:)];
end

%% naming a parameter file.

% Ask the user for a filename
current_dir = pwd; % Remember where we started
cd(shchpath)       % Move to shape change directory


% paramfilename=[current_shch(7:end) '.mat'];
% [~,paramfilenamebare,ext] = fileparts(paramfilename);

cd(current_dir)    % Go back to original directory
%

sysf_func = str2func(current_system);
shch_func = str2func(current_shch);
effnames = {'X','Y','Theta'};
if (OptFamily)
    optimizerdisplaytext = 'StepOpt: [';
    optimizerfiletext = 'stepopt_';
else
    optimizerdisplaytext = 'Opt: [';
    optimizerfiletext = 'opt_';
end
paramfiledisplaytext = [optimizerdisplaytext,sysf_func('name'),'] [',shch_func('name'),'] [',costfunction,'] ' effnames{direction} 'eff'];
paramfiletext = hash([optimizerfiletext,current_system(6:end),'_',current_shch(7:end),'_',costfunction,'_' effnames{direction} 'eff'],'md5');
paramfiletext = [optimizerfiletext,paramfiletext];

%% Save the data to a parameters file

% If the user didn't hit cancel, save the data and create a shchf file that
% reads the data and interprets it as a gait.
%  if ~usercancel

%     save(fullfile(shchpath,paramfilename),'alpha1','alpha2','t')
%     save(fullfile(shchpath,strcat(paramfilenamebare,'_optimal.mat')),'alpha1','alpha2','t')

if (n_dim >= 2) && (n_dim <= 5)
    % The optimal gait from step-optimizer is a 4x1 cell array.
    if(OptFamily)
        for i = 1:n_dim
            % The string makes "alphai = cell(4,1);"
            eval(['alpha',num2str(i),'=cell(4,1);']);
        end        
        for i = 1:n_dim
            for j = 1:4
                % The string makes "alphai{j} = alpha_out{j}(:,i);"
                eval(['alpha',num2str(i),'{',num2str(j),...
                    '} = alpha_out{',num2str(j),'}(:,',num2str(i),');']);
            end
        end        
    else        
        for i = 1:n_dim
            % The string makes "alphai = alpha_out(:,i);"
            eval(['alpha',num2str(i),'=alpha_out(:,',num2str(i),');']);
        end        
    end

    % save(fullfile(shchpath,strcat(paramfiletext,'.mat')),'alpha1','alpha2','t')
    save_string = "save(fullfile(shchpath,strcat(paramfiletext,'.mat')),";
    for i = 1:n_dim
        save_string = strcat(save_string,"'alpha",num2str(i),"',");
    end
    save_string = strcat(save_string, "'t')");
    eval(save_string);
        
else
    error('Trying to make an optimal gait with an unsupported number of dimensions')
end
    
    % Create the file if it doesn't already exist; future work could be
    % more fine-grained about this (making sure that any patches are
    % applied vs not overwriting any hand-edits the user made) and allowing
    % user to enter a prettier string for the display name here.
    
    % ['shchf_' paramfilenamebare '.m']
    if ~exist(fullfile(shchpath,['shchf_' paramfiletext '.m']),'file')
        % [paramfilenamebare '_optimal'], [paramfilenamebare '_optimal']
        gait_gui_draw_make_shchf(paramfiletext,paramfiledisplaytext,n_dim,OptFamily)
        
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