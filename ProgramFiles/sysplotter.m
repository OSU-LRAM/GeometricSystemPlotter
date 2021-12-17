function varargout = sysplotter(varargin)
% SYSPLOTTER M-file for sysplotter.fig
%      SYSPLOTTER, by itself, creates a new SYSPLOTTER or raises the existing
%      singleton*.
%
%      H = SYSPLOTTER returns the handle to a new SYSPLOTTER or the handle to
%      the existing singleton*.
%
%      SYSPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYSPLOTTER.M with the given input arguments.
%
%      SYSPLOTTER('Property','Value',...) creates a new SYSPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sysplotter_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sysplotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sysplotter

% Last Modified by GUIDE v2.5 06-May-2021 16:25:08

    addpath('./Utilities')

	%path to gui functions
	addpath(genpath(GetFullPath('sysplotter_gui_fcns')),...
            genpath(GetFullPath('sys_calcpath_fcns')),...
            genpath(GetFullPath('sys_calcsystem_fcns')),...
            genpath(GetFullPath('sys_draw_fcns')), ...
            genpath(GetFullPath('sys_update_fcns')),...
            genpath(GetFullPath('sysplotter_config_fcns')),...
            genpath(GetFullPath('Utilities')),...
            genpath(GetFullPath('Animation')),...
            genpath(GetFullPath('Physics')),...
            genpath(GetFullPath('Geometry')),...
            genpath(GetFullPath('GaitOptimization')),...
            genpath(GetFullPath('CoordinateOptimization')));

	%%%
	% Ensure that system files are properly accessible

	% Check for existence of sysplotter_config.mat data file
	configfile = fullfile(pwd,'sysplotter_config.mat');
	if exist(configfile,'file')

		% Load the configuration file
		load(configfile);

		% Verify that the configuration file points to a valid location
		v = verify_configdir(inputpath);

	end

	if ~exist(configfile,'file') || ~v % Automatically run the configuration 
		%code if config file is not present and valid

		%run the configuration function
		uiwait(sysplotter_config);

		% Recheck the existence of the sysplotter_config.mat file
		if exist(configfile,'file')

			% Load the configuration file
			load(configfile);

		else
			warning('Configuration file not found or generated')
			return
		end
	end




	% Verify that the target directory has the necessary subdirectories
	v = verify_configdir(inputpath);

	% Only continue if the location contents are valid
	if v
		% Add the paths to the system directories
		addpath(fullfile(inputpath, '/Systems'),fullfile(inputpath, '/Shape_Changes'));
        
        % Assign a variable in the base workspace that points to the user
        % files
        assignin('base','sysplotter_inputpath',inputpath);
        
	else
		warning('Configuration file points to a directory without the necessary subdirectories')
		return
	end



	% Ensure that the directory for storing processed data exists
	if ~exist(datapath,'dir')
		mkdir(datapath);
	end
	addpath(datapath)

	% Begin initialization code - DO NOT EDIT
	gui_Singleton = 1;
	gui_State = struct('gui_Name',       mfilename, ...
		'gui_Singleton',  gui_Singleton, ...
		'gui_OpeningFcn', @sysplotter_OpeningFcn, ...
		'gui_OutputFcn',  @sysplotter_OutputFcn, ...
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


% --- Executes just before sysplotter is made visible.
function sysplotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sysplotter (see VARARGIN)

% Choose default command line output for sysplotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Disable the plotting button unless a reasonable value has been set in
% the system and shape change menus
plotpublishpushbutton_enable_menu(handles)

%Get the configuration file, and extract the Colorpath
configfile = './sysplotter_config';
load(configfile,'Colorset');

% Create an empty progress bar in the progress bar panel
waitbar2a(0,handles.progresspanel,'waitbartext','Waiting for input',...
    'barcolor',Colorset.spot);

% Update the handles with the desired properties
load(configfile,'propertyfilepath');
GUI_property_loader(handles,propertyfilepath)
	
end


% --- Outputs from this function are returned to the command line.
function varargout = sysplotter_OutputFcn(hObject, eventdata, handles)
	% varargout  cell array for returning output args (see VARARGOUT);
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Get default command line output from handles structure
	varargout{1} = handles.output;

end
