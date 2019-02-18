function [] = sysf_config_gui_initvals(hObject,varargin)
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
handles.curr_system.String = varargin{1}{1};

guidata(hObject,handles); %Save data back to the hObject figure

end

