function varargout = osciplot_gui_v2(varargin)
%OSCIPLOT_GUI_V2 MATLAB code file for osciplot_gui_v2.fig
%      OSCIPLOT_GUI_V2, by itself, creates a new OSCIPLOT_GUI_V2 or raises the existing
%      singleton*.
%       
%       osciplot_gui_v2 expects this input:
%       ('experiment',recording number,channel,condition,[repetitions]).
%       e.g. osciplot_gui_v2('G17512',6,1,1,[1:10])
%
%       Make shure the dataset of the experiment is complete; a stimulus
%       and the corresponding trace (anadata)
% See also: GUIDE, GUIDATA, GUIHANDLES
% Code by Jan Everaert 2018
% Edit the above text to modify the response to help osciplot_gui_v2

% Last Modified by GUIDE v2.5 30-Nov-2018 19:54:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @osciplot_gui_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @osciplot_gui_v2_OutputFcn, ...
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


% --- Executes just before osciplot_gui_v2 is made visible.
function osciplot_gui_v2_OpeningFcn(hObject, eventdata, handles, varargin)
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
for ii = 1:length(repetition)
    handles.trace_original =  [handles.trace_original anadata(handles.ds,channel,cond,repetition(ii))'];
end
handles.trace_original = handles.trace_original';
handles.trace = handles.trace_original;
handles.stim = [];
% Temp solution for data sets with for the stimulus a cell (don't know why)
% Stimulus is constructed out of stim information present in the dataset
% (added by Jan 2018)
% extract both channels added by Marta, Oct 2018
if iscell(handles.ds.Stim.Waveform(cond).Samples)
    for channel = 1:numel(handles.ds.Stim.Waveform(cond,:)) %for binaural stimulus
        stim_construct = [];
        for ii = 1:length(handles.ds.Stim.Waveform(cond,channel).Samples)
            stim_part = handles.ds.Stim.Waveform(cond,channel).Samples{ii};
            stim_part_rep = handles.ds.Stim.Waveform(cond,channel).Nrep(ii);
            stim_construct = [stim_construct;repmat(stim_part,stim_part_rep,1)];
            stim =  repmat(stim_construct,length(repetition),1);
        end
        handles.stim(:,channel) = stim;
    end
else
    % Seems to be non optimal (Marta)
    for ii = 1:length(repetition)
        handles.stim =  [handles.stim handles.ds.Stim.Waveform(cond).Samples'];
    end
    
end
% handles stim is (2,X) array (if 2 channels)
handles.stim = handles.stim';
wf = handles.ds.Stim.Waveform(cond,:);
handles.DAchan = [wf(:).DAchan]';
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
legend(handles.DAchan)

% Initialization of handles user input fields
handles = init_gui(hObject,handles);

% Choose default command line output for osciplot_gui_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes osciplot_gui_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = osciplot_gui_v2_OutputFcn(hObject, eventdata, handles)
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
% hObject    handle to edit_tFilterlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.filter.lowerthreshold = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);


% Hints: get(hObject,'String') returns contents of edit_tFilterlow as text
%        str2double(get(hObject,'String')) returns contents of edit_tFilterlow as a double


% --- Executes during object creation, after setting all properties.
function filter_LL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tFilterlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_UL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tFilterup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.filter.upperthreshold = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit_tFilterup as text
%        str2double(get(hObject,'String')) returns contents of edit_tFilterup as a double


% --- Executes during object creation, after setting all properties.
function filter_UL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tFilterup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_extFilter.
function filter_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_extFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
plot_demand = 1; %$ plot wanted
trace_filtered = filter_and_plot(hObject,handles,plot_demand);
handles.trace = trace_filtered;

% update the auto Y limits
handles.y_limits_trace.auto.min = min(handles.trace);
handles.y_limits_trace.auto.max = max(handles.trace);

guidata(hObject, handles);


% --- Executes on selection change in popup_tTrigger.
function popup_tTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to popup_tTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.peak.channel = get(hObject,'Value');
guidata(hObject,handles);


% Hints: contents = cellstr(get(hObject,'String')) returns popup_tTrigger contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_tTrigger


% --- Executes during object creation, after setting all properties.
function popup_tTrigger_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_tTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_extTrigger.
function peak_detection_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_extTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
plot_demand = 1;
[peaks, stim_int, t_int] = peak_finder_and_plot(hObject,handles,plot_demand);
handles.peak.positions = peaks;
handles.stim_int = stim_int;
handles.t_int = t_int;
guidata(hObject, handles);



function peak_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.peak.threshold = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit_tTrigger as text
%        str2double(get(hObject,'String')) returns contents of edit_tTrigger as a double


% --- Executes during object creation, after setting all properties.
function peak_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function window_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tWsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tWsize as text
%        str2double(get(hObject,'String')) returns contents of edit_tWsize as a double

handles.window_size.time = str2double(get(hObject,'String'));
handles.window_size.samples = round(handles.window_size.time*handles.Fs);

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function window_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tWsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_extWsize.
function window_size_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_extWsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
window_size_and_plot(hObject,handles);


function edit_tYlim_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tYlim_min (see GCBO)
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
function edit_tYlim_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tYlim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_tYlim_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tYlim_min (see GCBO)
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
function edit_tYlim_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tYlim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in radiobutton_tYlim_auto.
function radiobutton_tYlim_auto_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_tYlim_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.y_limits_trace.auto_on = 1;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in radiobutton_tYlim_man.
function radiobutton_tYlim_man_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_tYlim_man (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.y_limits_trace.auto_on = 0;
% Update handles structure
guidata(hObject, handles);


function display_speed_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mSpeed as text
%        str2double(get(hObject,'String')) returns contents of edit_mSpeed as a double
handles.display_speed_ratio = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function display_speed_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function window_step_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mWstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mWstep as text
%        str2double(get(hObject,'String')) returns contents of edit_mWstep as a double
handles.window_step_size.time = str2double(get(hObject,'String'));
handles.window_step_size.samples = round(handles.window_step_size.time*handles.Fs);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function window_step_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mWstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minimal_display_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mMinDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mMinDisp as text
%        str2double(get(hObject,'String')) returns contents of edit_mMinDisp as a double
handles.minimal_display_time.time = str2double(get(hObject,'String'));
handles.minimal_display_time.samples = round(handles.minimal_display_time.time*handles.Fs);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function minimal_display_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mMinDisp (see GCBO)
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
% hObject    handle to edit_mName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mName as text
%        str2double(get(hObject,'String')) returns contents of edit_mName as a double
handles.file_specs.name = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function file_location_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mFolder as text
%        str2double(get(hObject,'String')) returns contents of edit_mFolder as a double
handles.file_specs.location = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function file_location_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_mfolder.
function pushbutton_mfolder_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = fileparts(fileparts(mfilename('fullpath')));
mfolder = uigetdir(path,'Select folder');
if all(mfolder)
    handles.file_specs.location = mfolder;
    h_mFolder = findobj('Tag','edit_mFolder');
    h_mFolder.String = mfolder;
    handles.file_specs.location = mfolder; 
    guidata(hObject, handles);
end

% --- Executes on button press in pushbutton20.
function save_movie_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = 'Yes';
if isfile(strcat(handles.file_specs.location,'\',handles.file_specs.name,'.avi'))
    answer = questdlg('A file with the same name already exists. Do you want to overwrite it?', ...
	'Movie file name', ...
	'Yes','No','No');
end
if strcmp(answer, 'Yes')
    save_movie(hObject,handles);
end

   
% --- Executes when selected object is changed in uibuttongroup13.
function uibuttongroup13_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup13 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

		% Get the values of the radio buttons in this group.
		handles.source_audio.stim = get(handles.radiobutton_mAudstim, 'Value');
		handles.source_audio.trace = get(handles.radiobutton_mAudtrace, 'Value');

% Update handles structure
guidata(hObject, handles);



function edit_sYlim_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sYlim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sYlim_max as text
%        str2double(get(hObject,'String')) returns contents of edit_sYlim_max as a double


% --- Executes during object creation, after setting all properties.
function edit_sYlim_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sYlim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sYlim_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sYlim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sYlim_min as text
%        str2double(get(hObject,'String')) returns contents of edit_sYlim_min as a double


% --- Executes during object creation, after setting all properties.
function edit_sYlim_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sYlim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sTrigger as text
%        str2double(get(hObject,'String')) returns contents of edit_sTrigger as a double


% --- Executes during object creation, after setting all properties.
function edit_sTrigger_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_sTrigger.
function pushbutton_sTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popup_sTrigger.
function popup_sTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_sTrigger contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sTrigger


% --- Executes during object creation, after setting all properties.
function popup_sTrigger_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_sWsize.
function pushbutton_sWsize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sWsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
window_size_and_plot(hObject,handles);



function edit_sWsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sWsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sWsize as text
%        str2double(get(hObject,'String')) returns contents of edit_sWsize as a double
handles.window_size_s.time = str2double(get(hObject,'String'));
handles.window_size_s.samples = round(handles.window_size.time*handles.Fs);



% --- Executes during object creation, after setting all properties.
function edit_sWsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sWsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
