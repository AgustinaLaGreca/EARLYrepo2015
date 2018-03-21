function [] = movie_preview(hObject,handles)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Asserts preview

% The output frame rate should not exeed 1000fps (for the preview)
if handles.display_speed_ratio*handles.Fs/handles.window_step_size.samples > 1000
    msgbox(strcat('The output frame rate should not exeed 1000fps, please increase the window step size or decrease the display speed ratio.'));

    return
end


%% Preview of the movie and sound
% -> Selection of the window_step_size and display_speed_ratio
% fast visualisation of the video
handles.display_freq = handles.display_speed_ratio*handles.Fs/handles.window_step_size.samples;
FV = handles.frames.video;

figure()
movie(FV,5,handles.display_freq/3)

% Update handles structure
guidata(hObject, handles);
end

