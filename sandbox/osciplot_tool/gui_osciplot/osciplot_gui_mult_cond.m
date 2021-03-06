function varargout = osciplot_gui_mult_cond(varargin)
%OSCIPLOT_GUI MATLAB code file for osciplot_gui.fig
%      OSCIPLOT_GUI is a user interface that can be used as a tool for
%      creation of moving images of traces. The aim is to simulate (at
%      least to a certain extend) the image that can be appreciated on an
%      oscilloscope. This interface calls several matlab scibts.
%       
%       osciplot_gui expects this input:
%       ('experiment',recording number,channel,[conditions],repetition).
%       e.g. osciplot_gui_mult_cond('G17512',6,1,[1:3],1)
%
%       Make shure the dataset of the experiment is complete; a stimulus
%       and the corresponding trace (anadata)
% See also: GUIDE, GUIDATA, GUIHANDLES
% Code by Jan Everaert 2018
% Edit the above text to modify the response to help osciplot_gui

% Last Modified by GUIDE v2.5 15-May-2018 12:32:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @osciplot_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @osciplot_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before osciplot_gui is made visible.
function osciplot_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% experiment = 'G17512';
% recording number = 6;
% channel = 1;
% cond = 1;
% repetition = 1;

experiment = varargin{1};
recording_number = varargin{2};
channel = varargin{3};
cond = varargin{4};
repetition = varargin{5};

handles.experiment = experiment;
handles.recording_number = recording_number;
handles.channel = channel;
handles.cond = cond;
handles.repetition = repetition;

handles.ds=read(dataset,experiment,recording_number);
handles.trace_original = [];
for ii = 1:length(cond)
    handles.trace_original =  [handles.trace_original anadata(handles.ds,channel,cond(ii),repetition)'];
end
handles.trace_original = handles.trace_original';
handles.trace = handles.trace_original;
handles.stim = [];
% Solution for data sets with for the stimulus a cell (don't know why)
% Stimulus is constructed out of stim information present in the dataset
% (added by Jan 2018)
if iscell(handles.ds.Stim.Waveform(cond(1)).Samples)
for kk = 1:length(cond)
    
    stim_construct = [];
    for ii = 1:length(handles.ds.Stim.Waveform(cond).Samples)
    stim_part = handles.ds.Stim.Waveform(cond).Samples{ii};
    stim_part_rep = handles.ds.Stim.Waveform(cond).Nrep(ii);
    stim_construct = [stim_construct;repmat(stim_part,stim_part_rep,1)];
    end
    
    handles.stim =  [handles.stim stim_construct];
end    
else
    
    for kk = 1:length(cond)
        handles.stim =  [handles.stim handles.ds.Stim.Waveform(cond(kk)).Samples'];
    end
end 


handles.stim = handles.stim';
handles.Fs = handles.ds.Fsam;
handles.t = (0:length(handles.trace)-1)./handles.Fs;
handles.nb_of_samples = length(handles.trace_original);



% initial plot
axes(handles.axes1)
plot(handles.t,handles.trace)
axis tight
axes(handles.axes2)
plot(handles.t,handles.stim)
axis tight

% initialization
handles = init_gui(hObject,handles);

% Choose default command line output for osciplot_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes osciplot_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = osciplot_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


function filter_LL_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.filter.lowerthreshold = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);


% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function filter_LL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_UL_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.filter.upperthreshold = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function filter_UL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton9.
function filter_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
plot_demand = 1; %$ plot wanted
[trace_filtered] = filter_and_plot(hObject,handles,plot_demand);
handles.trace = trace_filtered;

% update the auto Y limits
handles.y_limits_trace.auto.min = min(handles.trace);
handles.y_limits_trace.auto.max = max(handles.trace);

guidata(hObject, handles);



% --- Executes on button press in pushbutton14.
function peak_detection_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
plot_demand = 1;
[peaks] = peak_finder_and_plot(hObject,handles,plot_demand);
handles.peak.positions = peaks;
guidata(hObject, handles);



function peak_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.peak.threshold = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function peak_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function window_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double

handles.window_size.time = str2double(get(hObject,'String'));
handles.window_size.samples = round(handles.window_size.time*handles.Fs);

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function window_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton15.
function window_size_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
window_size_and_plot(hObject,handles);


function display_speed_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double
handles.display_speed_ratio = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function display_speed_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function window_step_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double
handles.window_step_size.time = str2double(get(hObject,'String'));
handles.window_step_size.samples = round(handles.window_step_size.time*handles.Fs);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function window_step_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minimal_display_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double
handles.minimal_display_time.time = str2double(get(hObject,'String'));
handles.minimal_display_time.samples = round(handles.minimal_display_time.time*handles.Fs);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function minimal_display_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton17.
function execute_movie_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[video_frames,audio_frames_stim,audio_frames_trace] = get_all_frames(hObject,handles);
handles.frames.video = video_frames;
handles.frames.audio.stim = audio_frames_stim;
handles.frames.audio.trace = audio_frames_trace;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton18.
function preview_movie_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
movie_preview(hObject,handles);


function file_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double
handles.file_specs.name = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function file_location_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double
handles.file_specs.location = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function file_location_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton20.
function save_movie_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_movie(hObject,handles);

   
% --- Executes when selected object is changed in uibuttongroup13.
function uibuttongroup13_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup13 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

		% Get the values of the radio buttons in this group.
		handles.source_audio.stim = get(handles.radiobutton_stim, 'Value');
		handles.source_audio.trace = get(handles.radiobutton_trace, 'Value');

% Update handles structure
guidata(hObject, handles);




function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.y_limits_trace.manual.min = get(hObject,'String');
if iscell(handles.y_limits_trace.manual.min)
    handles.y_limits_trace.manual.min = cellfun(@str2num,handles.y_limits_trace.manual.min);
elseif ischar(handles.y_limits_trace.manual.min)
    handles.y_limits_trace.manual.min = str2double(handles.y_limits_trace.manual.min);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.y_limits_trace.manual.max = get(hObject,'String');
if iscell(handles.y_limits_trace.manual.max)
    handles.y_limits_trace.manual.max = cellfun(@str2num,handles.y_limits_trace.manual.max);
elseif ischar(handles.y_limits_trace.manual.max)
    handles.y_limits_trace.manual.max = str2double(handles.y_limits_trace.manual.max);
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.y_limits_trace.auto_on = 1;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.y_limits_trace.auto_on = 0;
% Update handles structure
guidata(hObject, handles);




