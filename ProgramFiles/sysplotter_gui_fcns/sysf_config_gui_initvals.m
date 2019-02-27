function [] = sysf_config_gui_initvals(hObject)
% sysf_config_gui_initvals:
% Description: This function will read all of the values from the mat config file
% if it exists and will then populate the gui for that specific file
% Inputs:
% hObject: This is the object for the whole figure to exchange guidata
% varargin: This is the input cell array with elements....
%   1: curr_system: name of the current system
%   2: matFilePath: Full filename and path to the saved mat file

handles = guidata(hObject); % Pull handles struct from the figure object

%Initialize the title
handles.curr_system.String = handles.varargin{1};

%Get the path to the saved config file
matFilePath = handles.varargin{2};

%Get the system structure
s = load(matFilePath,'s');
s = s.s;

%%%%%
% Set all of the text inputs
systemfields = {{'geometry','linklengths'},{'geometry','length'},...
    {'physics','drag_coefficient'},{'physics','drag_ratio'}, ...
    {'grid_range'},{'tic_locs','x'},{'tic_locs','y'},...
    {'visual','grid'},{'singularity'},{'density','vector'},...
    {'density','scalar'},{'density','eval'},...
    {'density','metric_eval'},{'density','finite_element'},{'geometry','baseframe'}};
n = 3;%Precision to show in text
mat2strn = @(x) mat2str(x,n);
num2strn = @(x) num2str(x,n);
valuefunc = {mat2strn,num2strn,num2strn,num2strn,mat2strn,mat2strn,mat2strn,...
    @(x) mat2str(x{1}(1:3)),num2strn,mat2strn,mat2strn,mat2strn,mat2strn,...
    num2strn,@(x) isnumbertostring(x,n)};


%Create all of the guielements
guielements = cell(size(systemfields));

for k = 1:size(systemfields,2)
    guielements{k} = handles.(strjoin(systemfields{k},'_'));
end

% Loop through all of the array fields
for k = 1:size(systemfields,2)
   if isfieldnested(s,systemfields{k})
        guielements{k}.String = valuefunc{k}(getfield(s,systemfields{k}{:})); 
   else
        guielements{k}.String = 'None';
   end
end

%Save the values to be used outside of the init function
handles.text_inputs.systemfields = systemfields;
handles.text_inputs.guielements = guielements;



%%%%%%%%
%Populate all of the drop down fields
dropdownfields = {{'geometry','type'},{'geometry','baseframe'}};
dropdownoptions = {{'None';'n-link chain';'general curvature';'curvature basis';'curvature bases'},...
    {'None';'Use Number';'center';'com-mean';'tail';'centered';'midpoint tangent';'head';'head-tip'}};

%Create all of the guielements
dropdownguielements = cell(size(systemfields));

for k = 1:size(dropdownfields,2)
    dropdownguielements{k} = handles.(strjoin([dropdownfields{k},{'popup'}],'_'));
end

% Loop through all of the array fields
for k = 1:size(dropdownfields,2)
    dropdownguielements{k}.String = dropdownoptions{k};
    %Find where the loaded value is in popup menu
    if isfieldnested(s,dropdownfields{k})
        position = find(strcmp(dropdownoptions{k}, getfield(s,dropdownfields{k}{:})));
        if isempty(position)
            if isnumeric(getfield(s,dropdownfields{k}{:}))
                dropdownguielements{k}.Value = 2;
            else
                dropdownguielements{k}.Value = 1;
            end
        else
            dropdownguielements{k}.Value = position;
        end
   end
end


%Save the values to be used outside of the init function
handles.text_inputs.dropdownfields = dropdownfields;
handles.text_inputs.dropdownguielements = dropdownguielements;
handles.text_inputs.dropdownoptions = dropdownoptions;

%%%%%%% Need to figure out a better way to deal with these variables
%Deal with link_number and constraint_direction
link_numberfields = {'geometry','constraint_list','link_number'};

if isfieldnested(s,link_numberfields)
    handles.geometry_link_number.String = mat2str(vertcat(s.geometry.constraint_list.link_number));
else
    handles.geometry_link_number.String = 'None';
end

link_numberfields = {'geometry','constraint_list','constraint_direction'};

if isfieldnested(s,link_numberfields)
    handles.geometry_constraint_direction.Data = vertcat(s.geometry.constraint_list.constraint_direction);
else
    handles.geometry_constraint_direction.Data = -1;
end

%%%%%%%%
guidata(hObject,handles); %Save data back to the hObject figure

end

