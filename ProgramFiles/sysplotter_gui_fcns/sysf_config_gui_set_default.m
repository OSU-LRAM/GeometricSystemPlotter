function [] = sysf_config_gui_set_default(hObject)
%sysf_config_gui_set_default: 
% Description: This function will be called when the set default button is
% pushed
% Input:
% hObject: Whole figure object
% varargin: The variables that came in as arguments to the original
% function
handles = guidata(hObject);

sysfun = str2func(handles.varargin{1});
sysfun('reset');
sysf_config_gui_initvals(hObject)

end

