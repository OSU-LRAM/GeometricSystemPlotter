function varargout = sysf_config_gui(varargin)
% SYSF_CONFIG_GUI MATLAB code for sysf_config_gui.fig
%      SYSF_CONFIG_GUI, by itself, creates a new SYSF_CONFIG_GUI or raises the existing
%      singleton*.
%
%      H = SYSF_CONFIG_GUI returns the handle to a new SYSF_CONFIG_GUI or the handle to
%      the existing singleton*.
%
%      SYSF_CONFIG_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYSF_CONFIG_GUI.M with the given input arguments.
%
%      SYSF_CONFIG_GUI('Property','Value',...) creates a new SYSF_CONFIG_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sysf_config_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sysf_config_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sysf_config_gui

% Last Modified by GUIDE v2.5 26-Feb-2019 22:45:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sysf_config_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @sysf_config_gui_OutputFcn, ...
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


% --- Executes just before sysf_config_gui is made visible.
function sysf_config_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sysf_config_gui (see VARARGIN)

% Choose default command line output for sysf_config_gui
handles.output = hObject;

%Save the input arguments to be used later
handles.varargin = varargin;

% Update handles structure
guidata(hObject, handles);

%Initialize all of the gui elements
sysf_config_gui_initvals(hObject)

handles.edit_default_button.Callback = @(src, evnt) edit([varargin{1} '.m']);

handles.restore_defaults_button.Callback = @(src, evnt) sysf_config_gui_set_default(hObject);

handles.save_button.Callback = @(src, evnt) sysf_config_gui_savevals(hObject);



% UIWAIT makes sysf_config_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sysf_config_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in geometry_type_popup.
function geometry_type_popup_Callback(hObject, eventdata, handles)
% hObject    handle to geometry_type_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns geometry_type_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from geometry_type_popup


% --- Executes during object creation, after setting all properties.
function geometry_type_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometry_type_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in geo_type_q.
function geo_type_q_Callback(hObject, eventdata, handles)
% hObject    handle to geo_type_q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in geometry_baseframe_popup.
function geometry_baseframe_popup_Callback(hObject, eventdata, handles)
% hObject    handle to geometry_baseframe_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns geometry_baseframe_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from geometry_baseframe_popup


% --- Executes during object creation, after setting all properties.
function geometry_baseframe_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometry_baseframe_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function geometry_linklengths_Callback(hObject, eventdata, handles)
% hObject    handle to geometry_linklengths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of geometry_linklengths as text
%        str2double(get(hObject,'String')) returns contents of geometry_linklengths as a double


% --- Executes during object creation, after setting all properties.
function geometry_linklengths_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometry_linklengths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function geometry_length_Callback(hObject, eventdata, handles)
% hObject    handle to geometry_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of geometry_length as text
%        str2double(get(hObject,'String')) returns contents of geometry_length as a double


% --- Executes during object creation, after setting all properties.
function geometry_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometry_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in edit_default_button.
function edit_default_button_Callback(hObject, eventdata, handles)
% hObject    handle to edit_default_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton51.
function pushbutton51_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tic_locs_y_Callback(hObject, eventdata, handles)
% hObject    handle to tic_locs_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tic_locs_y as text
%        str2double(get(hObject,'String')) returns contents of tic_locs_y as a double


% --- Executes during object creation, after setting all properties.
function tic_locs_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tic_locs_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tic_locs_x_Callback(hObject, eventdata, handles)
% hObject    handle to tic_locs_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tic_locs_x as text
%        str2double(get(hObject,'String')) returns contents of tic_locs_x as a double


% --- Executes during object creation, after setting all properties.
function tic_locs_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tic_locs_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function singularity_Callback(hObject, eventdata, handles)
% hObject    handle to singularity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of singularity as text
%        str2double(get(hObject,'String')) returns contents of singularity as a double


% --- Executes during object creation, after setting all properties.
function singularity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singularity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton47.
function pushbutton47_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function grid_range_Callback(hObject, eventdata, handles)
% hObject    handle to grid_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_range as text
%        str2double(get(hObject,'String')) returns contents of grid_range as a double


% --- Executes during object creation, after setting all properties.
function grid_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function visual_grid_Callback(hObject, eventdata, handles)
% hObject    handle to visual_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of visual_grid as text
%        str2double(get(hObject,'String')) returns contents of visual_grid as a double


% --- Executes during object creation, after setting all properties.
function visual_grid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visual_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_vector_Callback(hObject, eventdata, handles)
% hObject    handle to density_vector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_vector as text
%        str2double(get(hObject,'String')) returns contents of density_vector as a double


% --- Executes during object creation, after setting all properties.
function density_vector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_vector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function density_finite_element_Callback(hObject, eventdata, handles)
% hObject    handle to density_finite_element (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_finite_element as text
%        str2double(get(hObject,'String')) returns contents of density_finite_element as a double


% --- Executes during object creation, after setting all properties.
function density_finite_element_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_finite_element (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function density_scalar_Callback(hObject, eventdata, handles)
% hObject    handle to density_scalar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_scalar as text
%        str2double(get(hObject,'String')) returns contents of density_scalar as a double


% --- Executes during object creation, after setting all properties.
function density_scalar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_scalar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function density_eval_Callback(hObject, eventdata, handles)
% hObject    handle to density_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_eval as text
%        str2double(get(hObject,'String')) returns contents of density_eval as a double


% --- Executes during object creation, after setting all properties.
function density_eval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function density_metric_eval_Callback(hObject, eventdata, handles)
% hObject    handle to density_metric_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_metric_eval as text
%        str2double(get(hObject,'String')) returns contents of density_metric_eval as a double


% --- Executes during object creation, after setting all properties.
function density_metric_eval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_metric_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function physics_drag_ratio_label_Callback(hObject, eventdata, handles)
% hObject    handle to physics_drag_ratio_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of physics_drag_ratio_label as text
%        str2double(get(hObject,'String')) returns contents of physics_drag_ratio_label as a double


% --- Executes during object creation, after setting all properties.
function physics_drag_ratio_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to physics_drag_ratio_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function physics_drag_coefficient_Callback(hObject, eventdata, handles)
% hObject    handle to physics_drag_coefficient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of physics_drag_coefficient as text
%        str2double(get(hObject,'String')) returns contents of physics_drag_coefficient as a double


% --- Executes during object creation, after setting all properties.
function physics_drag_coefficient_CreateFcn(hObject, eventdata, handles)
% hObject    handle to physics_drag_coefficient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function geometry_link_number_Callback(hObject, eventdata, handles)
% hObject    handle to geometry_link_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of geometry_link_number as text
%        str2double(get(hObject,'String')) returns contents of geometry_link_number as a double


% --- Executes during object creation, after setting all properties.
function geometry_link_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometry_link_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in restore_defaults_button.
function restore_defaults_button_Callback(hObject, eventdata, handles)
% hObject    handle to restore_defaults_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function physics_drag_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to physics_drag_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of physics_drag_ratio as text
%        str2double(get(hObject,'String')) returns contents of physics_drag_ratio as a double


% --- Executes during object creation, after setting all properties.
function physics_drag_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to physics_drag_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function geometry_baseframe_Callback(hObject, eventdata, handles)
% hObject    handle to geometry_baseframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of geometry_baseframe as text
%        str2double(get(hObject,'String')) returns contents of geometry_baseframe as a double


% --- Executes during object creation, after setting all properties.
function geometry_baseframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometry_baseframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function visual_grid_num_cells_Callback(hObject, eventdata, handles)
% hObject    handle to visual_grid_num_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of visual_grid_num_cells as text
%        str2double(get(hObject,'String')) returns contents of visual_grid_num_cells as a double


% --- Executes during object creation, after setting all properties.
function visual_grid_num_cells_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visual_grid_num_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton55.
function pushbutton55_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
