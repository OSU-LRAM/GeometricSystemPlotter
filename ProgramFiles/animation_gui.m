function varargout = animation_gui(varargin)
% ANIMATION_GUI MATLAB code for animation_gui.fig
%      ANIMATION_GUI, by itself, creates a new ANIMATION_GUI or raises the existing
%      singleton*.
%
%      H = ANIMATION_GUI returns the handle to a new ANIMATION_GUI or the handle to
%      the existing singleton*.
%
%      ANIMATION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANIMATION_GUI.M with the given input arguments.
%
%      ANIMATION_GUI('Property','Value',...) creates a new ANIMATION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before animation_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to animation_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help animation_gui

% Last Modified by GUIDE v2.5 07-Sep-2018 12:53:48

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @animation_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @animation_gui_OutputFcn, ...
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


% --- Executes just before animation_gui is made visible.
function animation_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to animation_gui (see VARARGIN)

% Choose default command line output for animation_gui
handles.main_gui = varargin{1};
handles.output = hObject;



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes animation_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% % --- Executes on button press in animatebutton_gui2.
% function animatebutton_Callback(hObject, eventdata, handles)
% % hObject    handle to animatebutton_gui2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)



function number_input_Callback(hObject, eventdata, handles)
% hObject    handle to number_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_input as text
%        str2double(get(hObject,'String')) returns contents of number_input as a double


% --- Executes during object creation, after setting all properties.
function number_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function framrate_Callback(hObject, eventdata, handles)
% hObject    handle to framrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framrate as text
%        str2double(get(hObject,'String')) returns contents of framrate as a double


% --- Executes during object creation, after setting all properties.
function framrate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Outputs from this function are returned to the command line.
function varargout = animation_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function duration_Callback(hObject, eventdata, handles)
% hObject    handle to duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of duration as text
%        str2double(get(hObject,'String')) returns contents of duration as a double


% --- Executes during object creation, after setting all properties.
function duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_CCFx_movie.
function checkbox_CCFx_movie_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CCFx_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CCFx_movie


% --- Executes on button press in checkbox_CCFy_movie.
function checkbox_CCFy_movie_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CCFy_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CCFy_movie


% --- Executes on button press in checkbox_CCFtheta_movie.
function checkbox_CCFtheta_movie_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CCFtheta_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CCFtheta_movie


% --- Executes on button press in checkbox_system_movie.
function checkbox_system_movie_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_system_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_system_movie


% --- Executes on button press in checkbox_vfieldx_movie.
function checkbox_vfieldx_movie_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vfieldx_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vfieldx_movie


% --- Executes on button press in checkbox_vfieldy_movie.
function checkbox_vfieldy_movie_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vfieldy_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vfieldy_movie


% --- Executes on button press in checkbox_vfieldtheta_movie.
function checkbox_vfieldtheta_movie_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vfieldtheta_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vfieldtheta_movie
