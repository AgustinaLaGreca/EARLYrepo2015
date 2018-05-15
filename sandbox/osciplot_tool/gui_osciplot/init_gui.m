function [handles] = init_gui(hObject,handles)
%[handles] = init_gui(hObject,handles)
%   Initializes all the values for the different parameters on the gui.
%   Input: the signal parameters that are allready in handles (struct)
%   Output: An updated struct handles

%% filter init

% LL and UL init values
handles.filter.lowerthreshold = 300;
handles.filter.upperthreshold = 3000;

% filtering itself (calculation) % Not needed
% axes(handles.axes1)
% plot_demand = 0;
% [trace_filtered] = filter_and_plot(hObject,handles,plot_demand);
% handles.trace = trace_filtered;

%% Y limits plots
handles.y_limits_trace.auto_on = 1;
handles.y_limits_trace.auto.min = min(handles.trace);
handles.y_limits_trace.auto.max = max(handles.trace);
handles.y_limits_trace.manual.min = min(handles.trace);
handles.y_limits_trace.manual.max = max(handles.trace);

%% Peak detection init

% Treshold values peakfinder
handles.peak.threshold = 0.0075;

% Peakfinder excecution
axes(handles.axes1)
plot_demand = 0;
[peaks] = peak_finder_and_plot(hObject,handles,plot_demand);
handles.peak.positions = peaks;

%% Window size init

% window size init value
handles.window_size.time = 0.015;
handles.window_size.samples = round(handles.window_size.time*handles.Fs);

%% Movie parameters init

handles.display_speed_ratio = 1;

handles.window_step_size.time = 0.010;
handles.window_step_size.samples = round(handles.window_step_size.time*handles.Fs);

handles.minimal_display_time.time = 0.015;
handles.minimal_display_time.samples = round(handles.minimal_display_time.time*handles.Fs);

%% Audio selection init

handles.source_audio.stim = 0;
handles.source_audio.trace = 1;

%% Movie save file specs init

handles.file_specs.name = strcat('movie','_',num2str(handles.experiment),'_',...
    num2str(handles.recording_number),'_',...
    num2str(handles.channel),'_',num2str(handles.cond),'_',...
    num2str(handles.repetition));
handles.file_specs.location = 'C:\EARLYrepo2015\sandbox\osciplot_tool';


%% Text on the gui
% data info set
set(handles.text53, 'String', num2str(handles.experiment));
set(handles.text54, 'String', num2str(handles.recording_number));
set(handles.text55, 'String', num2str(handles.channel));
set(handles.text56, 'String', num2str(handles.cond));
set(handles.text57, 'String', num2str(handles.repetition));

% Filter
set(handles.edit9, 'String', handles.filter.lowerthreshold);
set(handles.edit10, 'String', handles.filter.upperthreshold);

% Y-limits
set(handles.edit37, 'String', round(handles.y_limits_trace.manual.min,4));
set(handles.edit38, 'String', round(handles.y_limits_trace.manual.max,4));

% Peak detection
set(handles.edit21, 'String', handles.peak.threshold);

% Window size
set(handles.edit23, 'String', handles.window_size.time);

% Movie param
set(handles.edit33, 'String', handles.display_speed_ratio);
set(handles.edit35, 'String', handles.window_step_size.time);
set(handles.edit36, 'String', handles.minimal_display_time.time);

% Movie save
set(handles.edit28, 'String', handles.file_specs.name);
set(handles.edit29, 'String', handles.file_specs.location);


