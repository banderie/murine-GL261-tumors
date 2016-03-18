function varargout = segment_gui(varargin)
% SEGMENT_GUI MATLAB code for segment_gui.fig
%      SEGMENT_GUI, by itself, creates a new SEGMENT_GUI or raises the existing
%      singleton*.
%
%      H = SEGMENT_GUI returns the handle to a new SEGMENT_GUI or the handle to
%      the existing singleton*.
%
%      SEGMENT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENT_GUI.M with the given input arguments.
%
%      SEGMENT_GUI('Property','Value',...) creates a new SEGMENT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segment_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segment_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segment_gui

% Last Modified by GUIDE v2.5 19-Oct-2014 22:38:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segment_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @segment_gui_OutputFcn, ...
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

% --- Executes just before segment_gui is made visible.
function segment_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segment_gui (see VARARGIN)

addpath /Users/banderies/Code/Brain-Tumors

% Initialize variables
handles.z_slice = 1;
handles.threshold = 0;
handles.process_type = 'segment';

% Choose default command line output for segment_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes segment_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = segment_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in data_list.
function data_list_Callback(hObject, eventdata, handles)
% hObject    handle to data_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data_list

val = get(hObject, 'Value');
str = get(hObject, 'String');
switch str{val}
    % S1_G2_M1
    case 'S1_G2_M1_scan_1'
        g = 2;
        m = 1;
        timept = 1;
    case 'S1_G2_M1_scan_2'
        g = 2;
        m = 1;
        timept = 2;
    case 'S1_G2_M1_scan_3'
        g = 2;
        m = 1;
        timept = 3;
    case 'S1_G2_M1_scan_4'
        g = 2;
        m = 1;
        timept = 4;
    case 'S1_G2_M1_scan_5'
        g = 2;
        m = 1;
        timept = 5;
        % S1_G2_M2
    case 'S1_G2_M2_scan_1'
        g = 2;
        m = 2;
        timept = 1;
    case 'S1_G2_M2_scan_2'
        g = 2;
        m = 2;
        timept = 2;
    case 'S1_G2_M2_scan_3'
        g = 2;
        m = 2;
        timept = 3;
    case 'S1_G2_M2_scan_4'
        g = 2;
        m = 2;
        timept = 4;
    case 'S1_G2_M2_scan_5'
        g = 2;
        m = 2;
        timept = 5;
        % S1_G3_M1
    case 'S1_G3_M1_scan_1'
        g = 3;
        m = 1;
        timept = 1;
    case 'S1_G3_M1_scan_2'
        g = 3;
        m = 1;
        timept = 2;
    case 'S1_G3_M1_scan_3'
        g = 3;
        m = 1;
        timept = 3;
    case 'S1_G3_M1_scan_4'
        g = 3;
        m = 1;
        timept = 4;
    case 'S1_G3_M1_scan_5'
        g = 3;
        m = 1;
        timept = 5;
        % S1_G3_M2
    case 'S1_G3_M2_scan_1'
        g = 3;
        m = 2;
        timept = 1;
    case 'S1_G3_M2_scan_2'
        g = 3;
        m = 2;
        timept = 2;
    case 'S1_G3_M2_scan_3'
        g = 3;
        m = 2;
        timept = 3;
    case 'S1_G3_M2_scan_4'
        g = 3;
        m = 2;
        timept = 4;
    case 'S1_G3_M2_scan_5'
        g = 3;
        m = 2;
        timept = 5;
        % S1_G3_M3
    case 'S1_G3_M3_scan_1'
        g = 3;
        m = 3;
        timept = 1;
    case 'S1_G3_M3_scan_2'
        g = 3;
        m = 3;
        timept = 2;
    case 'S1_G3_M3_scan_3'
        g = 3;
        m = 3;
        timept = 3;
    case 'S1_G3_M3_scan_4'
        g = 3;
        m = 3;
        timept = 4;
    case 'S1_G3_M3_scan_5'
        g = 3;
        m = 3;
        timept = 5;
end

res = 0.5;

handles.filename_brain = experiment_filename(g,m,res,timept,'brain');
handles.filename_ventricle = experiment_filename(g,m,res,timept,'ventricles');
handles.filename_tumor = experiment_filename(g,m,res,timept,'tumor all');

%disp(handles.filename_brain);
%disp(handles.filename_ventricle);
%disp(handles.filename_tumor);

process (hObject, eventdata, handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function data_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in z_slice_textbox.
function z_slice_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to z_slice_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns z_slice_textbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from z_slice_textbox

handles.z_slice = get(hObject, 'Value');

guidata(hObject, handles);

process (hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function z_slice_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_slice_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threshold_text_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_text as text
%        str2double(get(hObject,'String')) returns contents of threshold_text as a double


handles.threshold = get(hObject, 'String');
handles.threshold = str2num(handles.threshold);

% disp(handles.threshold)

guidata(hObject, handles);

process (hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function threshold_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function threshold_brain (hObject, eventdata, handles)

filenames.brain = handles.filename_brain;
filenames.tumor = handles.filename_tumor;
filenames.ventricle = handles.filename_ventricle;

[~,~,raw,~,~,~] = mask(filenames);

mod_raw = raw;
% min_density = min(min(min(raw)));
% max_density = max(max(max(raw)));

% threshold = 6600;
mod_raw(mod_raw < handles.threshold) = 0;

image(squeeze(raw(:,:,handles.z_slice)), 'Parent', handles.axes1);
colormap(bone(20000));

image(squeeze(mod_raw(:,:,handles.z_slice)), 'Parent', handles.axes2);
colormap(bone(20000));

function segment_tumor (hObject, eventdata, handles)

[raw,tumor,ventricle,x,y,z] = loadData(handles.filename_brain,handles.filename_tumor,handles.filename_ventricle);

mod_raw = raw;

for k = 1:z
    for j = 1:y
        for i = 1:x
            if ~isnan(ventricle(i,j,k))
                raw(i,j,k) = NaN;
            end
        end
    end
end

for k = 1:z
    for j = 1:y
        for i = 1:x
            if isnan(tumor(i,j,k))
                mod_raw(i,j,k) = NaN;
            end
        end
    end
end

image(squeeze(raw(:,:,handles.z_slice)), 'Parent', handles.axes1);
colormap(bone(20000));

image(squeeze(mod_raw(:,:,handles.z_slice)), 'Parent', handles.axes2);
colormap(bone(20000));


% --- Executes on button press in threshold_pushbutton.
function threshold_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.process_type = 'threshold';
guidata(hObject, handles);
process(hObject, eventdata, handles);

% --- Executes on button press in segment_pushbutton.
function segment_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to segment_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of segment_pushbutton

handles.process_type = 'segment';
guidata(hObject, handles);
process(hObject, eventdata, handles);

function process (hObject, eventdata, handles)

if strcmp(handles.process_type, 'threshold')
    threshold_brain (hObject, eventdata, handles);
elseif strcmp(handles.process_type,'segment')
    segment_tumor (hObject, eventdata, handles)
end
    