function [video_frames,audio_frames_stim,audio_frames_trace] = get_all_frames(hObject,handles)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Asserts get all frames

% Display speed ratio should be strictly positive
if handles.display_speed_ratio <= 0
    msgbox(strcat('The display speed ratio should be striclty positive, please change your input.'));
    
    video_frames = [];audio_frames_stim=[];audio_frames_trace=[];
    return
end

% Window step size should be larger than one time step (1/Fs) and smaller
% than half the signal size
if handles.window_step_size.time < 1/handles.Fs
    msgbox(strcat('The minimal window step size should be larger than one time step (=',num2str(1/handles.Fs),'s), please change your input.'));
    
    video_frames = [];audio_frames_stim=[];audio_frames_trace=[];
    return
end

if handles.window_step_size.samples > handles.nb_of_samples/2
    msgbox(strcat('The step size should be smaller than half the signal size (=',num2str(handles.nb_of_samples/2/handles.Fs),'s), please change your input.'));
    
    video_frames = [];audio_frames_stim=[];audio_frames_trace=[];
    return
end

% The minimal display time should be positive or zero
if handles.minimal_display_time.samples < 0
    msgbox(strcat('The minimal display time should be positive or zero, please change your input.'));
    
    video_frames = [];audio_frames_stim=[];audio_frames_trace=[];
    return
end

% The output frame rate should not exeed 1000fps (for the preview)
if handles.display_speed_ratio*handles.Fs/handles.window_step_size.samples > 1000
    msgbox(strcat('The output frame rate should not exeed 1000fps, please increase the window step size or decrease the display speed ratio.'));
    
    video_frames = [];audio_frames_stim=[];audio_frames_trace=[];
    return
end


%% get frames

stim = handles.stim;
trace = handles.trace;
peak_locs = handles.peak.positions;
window_step_size = handles.window_step_size.samples;
window_size = handles.window_size.samples;
nb_frames_fading = 15;
nb_samples = length(trace);
display_speed_ratio = handles.display_speed_ratio;
Fs = handles.Fs;
t_min_disp = handles.minimal_display_time.samples;

% fix the plot window (fixed axis)
max_stim = max(stim);
min_stim = min(stim);
max_trace = max(trace);
min_trace = min(trace);

alpha = 0.80; % Fading parameter

figure()
k = 1; % iteration parameter for the savings of the frames in each step
current_peaknr = 1;
if isempty(peak_locs) % when no peak is selected by the peakfinder
    current_peak_ind = 0;
else
    current_peak_ind = peak_locs(current_peaknr); 
end

for ii = 1:window_step_size:nb_samples-window_size
    

    frame_trace = trace(ii:ii+window_size-1);
    frame_stim = stim(ii:ii+window_size-1);
    
    % Save the corresponding audio frame
    FA_stim(:,k) = stim(ii:ii+window_step_size-1); %audio stim
    FA_trace(:,k) = trace(ii:ii+window_step_size-1); %audio trace

    
    if  current_peaknr < length(peak_locs)
        % go to the next peak, if there is one
        next_peaknr = current_peaknr+1;
        next_peak_ind = peak_locs(next_peaknr);

        if (next_peak_ind < ii+window_size)...
            && (ii > window_size/2) && (ii < nb_samples-window_size-window_size/2)
        
            current_peaknr = next_peaknr;
            current_peak_ind = peak_locs(current_peaknr);
        end
    end
    
    if (current_peak_ind < ii+window_size) && (ii-current_peak_ind)/(display_speed_ratio) < t_min_disp ...
            && (current_peak_ind > window_size/2) && (current_peak_ind < nb_samples-window_size-window_size/2)
        % peak is centered in the frame when it is in it
        new_frame_ind = ((current_peak_ind-floor((window_size)/2)):1:...
             (current_peak_ind+floor((window_size-1)/2)));
        frame_trace = trace(new_frame_ind);   
    end
    

    
    M(:,k) = frame_trace; % Save the frames for the fading
    N(:,k) = frame_stim; % not really used yet
   
    % Fade out of previous displays
    if k <= nb_frames_fading  
        subplot(211)
        for ll = 1:k-1
        plot(M(:,ll),'color',[0,0,0]+1-(alpha^(k-ll))^2,'LineWidth',1.2);hold on;
        end
    else 
        subplot(211)
        for ll = k-nb_frames_fading:k-1
        plot(M(:,ll),'color',[0,0,0]+1-(alpha^(k-ll))^2,'LineWidth',1.2);hold on;
        end 
    
    % Current display
    plot(M(:,k),'color',[0,0,0],'LineWidth',1.2);hold on %trace
    ylim([min_trace max_trace]);
    hold off;
    end
    
    subplot(212)
    plot(frame_stim); %stim
    ylim([min_stim max_stim]);
    
    % Save complete frame for video
    FV(k) = getframe(gcf);

    
    k=k+1;
end

video_frames = FV;
audio_frames_stim = FA_stim;
audio_frames_trace = FA_trace;
end

