function [] = save_movie(hObject,handles)
% [] = save_movie(hObject,handles)
%   Saves the constructed movie with the selected name
%   (handles.file_specs.name), in the selected location
%   (handles.file_specs.location). 
%   Make sure that the if you want to overwrite a movie (save in same place and with the same name as an
%   existing movie) that this movie is not opened by another program at the
%   same time

%% Actual video and audio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write audio and link with video frames and save as .avi

% kill WMP to be shure the .avi file can be overwritten
% If movie is played by another program, that should be killed too
% !taskkill -f -im wmplayer.exe

display_freq = handles.display_speed_ratio*handles.Fs/handles.window_step_size.samples;
FV = handles.frames.video;

if handles.source_audio.stim == 1
    FA = handles.frames.audio.stim;
else
    FA = handles.frames.audio.trace;
    FA = FA.*1000; %amplification of the signal needed
end
% construct movie and combine with audio and save it to a certain path
% if path is not specified it will be saved in the current folder (opened in matlab)

writerObj = vision.VideoFileWriter(strcat(handles.file_specs.location,'\',handles.file_specs.name,'.avi')...
    ,'AudioInputPort',true);
% writerObj = vision.VideoFileWriter('C:\EARLYrepo2015\sandbox\osciplot_tool\peaks2.avi','AudioInputPort',true);
writerObj.FrameRate = display_freq;
writerObj.VideoCompressor ='None (uncompressed)'; 
nb_of_frames = length(FV);
nb_of_replays_video = 1;

%%% Downsampling %%% If display freq is really high and the .avi file to large
% display_freq;
% downsample_rate = 3;
% writerObj.FrameRate = display_freq/downsample_rate;
% FV = downsample(FV,downsample_rate);
% FA = downsample(FA,downsample_rate);
% nb_of_frames = length(FV);



for ll = 1:nb_of_replays_video
    for i = 1 :1: nb_of_frames
        Frame=FV(i).cdata;
        %Frame = rgb2gray(Frame); %turn frames grayscale
        % add audio frames
        step(writerObj,Frame,FA(:,i));
    end
end
release(writerObj)

msgbox('Movie successfully saved')
end

