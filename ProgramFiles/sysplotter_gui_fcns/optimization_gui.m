function varargout = optimization_gui(varargin)
% OPTIMIZATION_GUI MATLAB code for optimization_gui.fig
%      OPTIMIZATION_GUI, by itself, creates a new OPTIMIZATION_GUI or raises the existing
%      singleton*.
%
%      H = OPTIMIZATION_GUI returns the handle to a new OPTIMIZATION_GUI or the handle to
%      the existing singleton*.
%
%      OPTIMIZATION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTIMIZATION_GUI.M with the given input arguments.
%
%      OPTIMIZATION_GUI('Property','Value',...) creates a new OPTIMIZATION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before optimization_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to optimization_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help optimization_gui

% Last Modified by GUIDE v2.5 22-Mar-2022 21:25:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @optimization_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @optimization_gui_OutputFcn, ...
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


% --- Executes just before optimization_gui is made visible.
function optimization_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to optimization_gui (see VARARGIN)

handles.main_gui = varargin{1};
% Choose default command line output for optimization_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes optimization_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = optimization_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in optimize_pushbutton.
function optimize_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to optimize_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the last plot pushbutton used
if isfield(handles.main_gui.figure1.UserData,'lastpushbutton')
    lastpushbutton = handles.main_gui.figure1.UserData.lastpushbutton;
else
    lastpushbutton = 'plotpushbutton1';
end

plot_info = plotpushbutton_Callback_optimize(findall(0,'tag',lastpushbutton), eventdata, handles.main_gui);

% Execute the gait_gui_optimize command
gait_gui_optimize(plot_info(1).axes(1),hObject, eventdata, handles, 0);
waitbar2a(1,handles.main_gui.progresspanel,'Finished Plotting')


% --- Executes on button press in stepoptimize_pushbutton.
function stepoptimize_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to stepoptimize_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the last plot pushbutton used
if isfield(handles.main_gui.figure1.UserData,'lastpushbutton')
    lastpushbutton = handles.main_gui.figure1.UserData.lastpushbutton;
else
    lastpushbutton = 'plotpushbutton1';
end

plot_info = plotpushbutton_Callback_optimize(findall(0,'tag',lastpushbutton), eventdata, handles.main_gui);

% Execute the gait_gui_optimize command
gait_gui_optimize(plot_info(1).axes(1),hObject, eventdata, handles, 1);
waitbar2a(1,handles.main_gui.progresspanel,'Finished Plotting')

% --- Executes on button press in thetabutton.
function thetabutton_Callback(hObject, eventdata, handles)
% hObject    handle to thetabutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of thetabutton


% --- Executes on button press in steeringbutton.
function steeringbutton_Callback(hObject, eventdata, handles)
% hObject    handle to steeringbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of steeringbutton


% --- Executes on button press in xbutton.
function xbutton_Callback(hObject, eventdata, handles)
% hObject    handle to xbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of xbutton


% --- Executes on button press in ybutton.
function ybutton_Callback(hObject, eventdata, handles)
% hObject    handle to ybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ybutton
