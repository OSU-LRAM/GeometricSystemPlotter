function varargout = sysplotter_config(varargin)
% SYSPLOTTER_CONFIG M-file for sysplotter_config.fig
%      SYSPLOTTER_CONFIG, by itself, creates a new SYSPLOTTER_CONFIG or raises the existing
%      singleton*.
%
%      H = SYSPLOTTER_CONFIG returns the handle to a new SYSPLOTTER_CONFIG or the handle to
%      the existing singleton*.
%
%      SYSPLOTTER_CONFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYSPLOTTER_CONFIG.M with the given input arguments.
%
%      SYSPLOTTER_CONFIG('Property','Value',...) creates a new SYSPLOTTER_CONFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sysplotter_config_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sysplotter_config_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sysplotter_config

% Last Modified by GUIDE v2.5 02-Dec-2015 17:46:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sysplotter_config_OpeningFcn, ...
                   'gui_OutputFcn',  @sysplotter_config_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

end

% --- Executes just before sysplotter_config is made visible.
function sysplotter_config_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sysplotter_config (see VARARGIN)

% Choose default command line output for sysplotter_config
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sysplotter_config wait for user response (see UIRESUME)
% uiwait(handles.sysplotter_config_gui);

% Add a path for subfunctions and utilities
addpath('sysplotter_config_fcns/');
addpath('Utilities');
    
	% Default values for HH and refpoint paths
    set(handles.inputpathconfig,'String',GetFullPath('../UserFiles/GenericUser'))
	set(handles.Colorconfig,'String','sys_draw_fcns/colorsets/color_Red.m')
	

    % Check for existence of sysplotter_config.mat data file
    configfile = './sysplotter_config.mat';
    if exist(configfile,'file')
        load(configfile);
        set(handles.inputpathconfig,'String',inputpath)
        set(handles.Colorconfig,'String',Colorpath)
    end
    
    % Verify that the target directory has the necessary subdirectories
    targetdir = get(handles.inputpathconfig,'String');
    v = verify_configdir(targetdir);

    % Only activate the OK button if the string is valid
    if v
        set(handles.OKbutton,'enable','on');
    else
        set(handles.OKbutton,'enable','off');
    end
    % Check version and set default property file
    if ispc
        propertyFile = 'sysplotter_config_fcns/propertyDataWindows.mat';
    elseif ismac
        propertyFile = 'sysplotter_config_fcns/propertyDataMacOS.mat';
    elseif isunix
        propertyFile = 'sysplotter_config_fcns/propertyDataLinux.mat';
    end
    set(handles.displayConfigFile,'String',propertyFile)

end


% --- Outputs from this function are returned to the command line.
function varargout = sysplotter_config_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end
