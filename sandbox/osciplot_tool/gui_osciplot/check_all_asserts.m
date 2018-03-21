function [] = check_all_asserts(hObject,handles);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% filter asserts

% LL and UL should be positive
if handles.filter.lowerthreshold <= 0 || handles.filter.upperthreshold <= 0
    msgbox('Please insert strictly positive values for the cuttof frequencies of the filter, please change your input.');
end

% LL should be smaller than UL
if handles.filter.lowerthreshold >= handles.filter.upperthreshold 
    msgbox('The lower CUF should be smaller than the upper CUF, please change your input.');
end


% UL should be smaller than the Nyquist frequency (=0.5*Fs)
if handles.filter.upperthreshold >= handles.Fs/2
    msgbox(strcat('The upper CUF should be smaller than the Nyquist freq (=',num2str(floor(handles.Fs/2)),'Hz), please change your input.'));
end


%% Peak detection asserts

% Threshold should be positive
if handles.peak.threshold < 0
    msgbox('Please insert a positive value for the threshold, please change your input.');
end


%% Window size

% Window size should be larger than 1 sample (time) = 1/Fs
if handles.window_size.samples < 1
    msgbox(strcat('The window size should be larger than one time step (=',num2str(1/handles.Fs),'s), please change your input.'));
end
% Window size should be smaller than the total signal (without margin) =
% length(trace)*Fs - window size
if handles.window_size.samples > length(handles.trace)/2
    msgbox(strcat('The window size should be smaller than half the total signal (=',num2str(length(handles.trace)/2/handles.Fs),'s), please change your input.'));                                                             
end

%% Display speed ratio

% DSR should be strictly positive
if handles.display_speed_ratio <= 0
    msgbox(strcat('The display speed ratio should be striclty positive, please change your input.'));
end

