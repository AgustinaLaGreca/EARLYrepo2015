function [handles] = init_gui(hObject,handles)
%[handles] = init_gui(hObject,handles)
%   Initializes all the values for the different user input parameters on
%   the gui (adds them to handles.
%   Input: the signal parameters that are allready in handles (struct)
%   Output: An updated struct handles

%% Get uicontrol objects
% TRACE -------------------------
% Filter
h_tFilterLow = findobj('Tag','edit_tFilterlow');
h_tFilterUp = findobj('Tag','edit_tFilterup');

% Trigger
hpop_tTrigger = findobj('Tag','popup_tTrigger');
hedit_tTrigger = findobj('Tag','edit_tTrigger');

% Window size
h_tWsize = findobj('Tag','edit_tWsize');

% Ylim
hbutton_tYlim = findobj('Tag','radiobutton_tYlim_auto');
h_tYlim_max = findobj('Tag','edit_tYlim_max');
h_tYlim_min = findobj('Tag','edit_tYlim_min');


% MOVIE -------------------------------
h_mSpeed = findobj('Tag','edit_mSpeed');
h_mWstep = findobj('Tag','edit_mWstep');
h_mDisp = findobj('Tag','edit_mMinDisp');

h_mAudTrace = findobj('Tag','radiobutton_mAudtrace');

h_mFilename = findobj('Tag','edit_mName');
h_mFolder = findobj('Tag','edit_mFolder');

%% filter init

handles.filter.lowerthreshold = str2double(get(h_tFilterLow,'String'));
handles.filter.upperthreshold = str2double(get(h_tFilterUp,'String'));

%% Y limits plots

handles.y_limits_trace.auto_on = (get(hbutton_tYlim,'Value'));
handles.y_limits_trace.auto.min = min(handles.trace);
handles.y_limits_trace.auto.max = max(handles.trace);
handles.y_limits_trace.manual.min = min(handles.trace);
handles.y_limits_trace.manual.max = max(handles.trace);

% Update values on GUI
h_tYlim_max.String = handles.y_limits_trace.manual.max;
h_tYlim_min.String = handles.y_limits_trace.manual.min;

%% Peak detection init
% Update popup menu with Experiment channels
str = ['Trace'; cellstr(handles.DAchan)];
hpop_tTrigger.String = str;

% Channel to used as trigger
handles.peak.channel = get(hpop_tTrigger,'Value');

% Threshold values peakfinder
handles.peak.threshold = str2double(get(hedit_tTrigger,'String'));

% Peakfinder excecution
axes(handles.axes1)
plotting = 0;
peaks = peak_finder_and_plot(hObject,handles,plotting);
handles.peak.positions = peaks;

%% Window size init

handles.window_size.time = str2double(get(h_tWsize,'String'));
handles.window_size.samples = round(handles.window_size.time*handles.Fs);

%% Movie parameters init

handles.display_speed_ratio = str2double(get(h_mSpeed,'String'));

handles.window_step_size.time = str2double(get(h_mWstep,'String'));
handles.window_step_size.samples = round(handles.window_step_size.time*handles.Fs);

handles.minimal_display_time.time = str2double(get(h_mDisp,'String'));
handles.minimal_display_time.samples = round(handles.minimal_display_time.time*handles.Fs);

%% Audio selection init

handles.source_audio.stim = ~(get(h_mAudTrace,'Value'));
handles.source_audio.trace = (get(h_mAudTrace,'Value'));

%% Movie save file specs init

handles.file_specs.name = strcat('movie','_',num2str(handles.experiment),'_',...
    num2str(handles.recording_number),'_',...
    num2str(handles.channel),'_',num2str(handles.cond),'_',...
    num2str(handles.repetition));
h_mFilename.String = handles.file_specs.name;
handles.file_specs.location = get(h_mFolder,'String');

%% Text on the gui
% data info set
set(handles.text53, 'String', num2str(handles.experiment));
set(handles.text54, 'String', num2str(handles.recording_number));
set(handles.text55, 'String', num2str(handles.channel));
set(handles.text56, 'String', num2str(handles.cond));
set(handles.text57, 'String', num2str(handles.repetition));

% Filter
%set(handles.edit9, 'String', handles.filter.lowerthreshold);
%set(handles.edit10, 'String', handles.filter.upperthreshold);

% Y-limits
%set(handles.edit_tYlim_min, 'String', round(handles.y_limits_trace.manual.min,4));
%set(handles.edit_tYlim_max, 'String', round(handles.y_limits_trace.manual.max,4));

% Peak detection
%set(handles.edit21, 'String', handles.peak.threshold);

% Window size
%set(handles.edit23, 'String', handles.window_size.time);

% Movie param
%set(handles.edit33, 'String', handles.display_speed_ratio);
%set(handles.edit35, 'String', handles.window_step_size.time);
%set(handles.edit36, 'String', handles.minimal_display_time.time);

% Movie save
%set(handles.edit28, 'String', handles.file_specs.name);
%set(handles.edit29, 'String', handles.file_specs.location);


