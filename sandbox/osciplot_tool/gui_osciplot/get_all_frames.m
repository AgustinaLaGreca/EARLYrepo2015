function [video_frames,audio_frames_stim,audio_frames_trace] = get_all_frames(hObject,handles)
% [video_frames,audio_frames_stim,audio_frames_trace] = get_all_frames(hObject,handles)
% returns the video frames and audio frames for the movie, based on all the
% selected parameters (stored in handles).

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

% For optimal performace, the window step size should be smaller than 2/3 of the
% window size

if handles.window_step_size.samples > handles.window_size.samples*2/3
    msgbox(strcat('For optimal visual performace, the step size should be smaller than 2/3  of the the window size (=',num2str(handles.window_size.samples*2/3/handles.Fs),'s), please change your input.'));
    
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
% If channeling in stim (peak.channel == 2 | == 3)
% use interpolated time and stimulus
if handles.peak.channel > 1
    time = handles.t_int;
    stim = handles.stim_int;
else
    time = handles.t;
    stim = handles.stim;
end
trace = handles.trace;
peak_locs = handles.peak.positions;
window_step_size = handles.window_step_size.samples;
window_size = handles.window_size.samples;
nb_frames_fading = 10;
nb_samples = length(trace);
display_speed_ratio = handles.display_speed_ratio;
Fs = handles.Fs;
t_min_disp = handles.minimal_display_time.samples;

% fix the plot window (fixed axis)
max_stim = max(stim(:)); 
min_stim = min(stim(:));
if handles.y_limits_trace.auto_on
    y_min_trace = handles.y_limits_trace.auto.min;
    y_max_trace = handles.y_limits_trace.auto.max;
else
    y_min_trace = handles.y_limits_trace.manual.min;
    y_max_trace = handles.y_limits_trace.manual.max;
end

alpha = 0.82; % Fading parameter

framesfig = figure('CloseRequestFcn',@close_showoff);
wb = waitbar(0,'0%','Name','Creating frames in background');
%k = 1; % iteration parameter for the savings of the frames in each step
%current_peaknr = 1;
if isempty(peak_locs) % when no peak is selected by the peakfinder (peakfinder_and_plot.m)
    current_peak_ind = 0;
else
    current_peak_ind = peak_locs(1); 
end

ind_frames = 1:window_step_size:nb_samples-window_size;
M = zeros(window_size,numel(ind_frames));
for k=1:numel(ind_frames) 
    ii = ind_frames(k);
    frame_trace = trace(ii:ii+window_size-1);
    frame_stim = stim(:,ii:ii+window_size-1);
    frame_time = time(ii:ii:ii+window_size-1);
    
    % Save the corresponding audio frame
    FA_stim(:,k) = stim(1,ii:ii+window_step_size-1); %audio stim - only one channel
    FA_trace(:,k) = trace(ii:ii+window_step_size-1); %audio trace
    frame_ind = ii:ii+window_step_size-1; % used for x-axis in the plot
    
    % get current peak nb (in current window second half, but not in next window)
    peaks_in_interval = peak_locs(peak_locs > ii & peak_locs < ii+window_step_size-1);
    if ~isempty(peaks_in_interval) && (ii > window_size/2) && (ii < nb_samples-window_size-window_size/2)
    current_peak_ind = peaks_in_interval(1);
    end
    
    % peak is centered in the frame when current frame contains a peak
    % first frames without peaks are not centered
    % if peak is not in current window (ie it was in previous window) and time between current window and previous spike is
    % shorter than t_min_disp, the displayed frame will continue centered
    % in previous spike. If it is larger, the frame is current window (no
    % centering done).
    %if handles.peak.channel == 1
    if (current_peak_ind < ii+window_size) && (ii-current_peak_ind)/(display_speed_ratio) < t_min_disp ...
            && (current_peak_ind > window_size/2) && (current_peak_ind < nb_samples-window_size-window_size/2)
        
        new_frame_ind = ((current_peak_ind-floor((window_size)/2)):1:...
             (current_peak_ind+floor((window_size-1)/2)));
        
        frame_trace = trace(new_frame_ind);
        frame_stim = stim(:,new_frame_ind);
        frame_time = time(new_frame_ind);
    %end
    end
        
    M(:,k) = frame_trace; % Save the frames for the fading
    %N(:,k) = frame_stim; % not really used yet
   
    % Update waitbar
    waitbar(k/numel(ind_frames),wb, [sprintf('%12.2f',k/numel(ind_frames)*100) '%']);
    % Fade out of previous displays
    set(0, 'CurrentFigure', framesfig);
    subplot(211)
    if k <= nb_frames_fading  
        for ll = 1:k-1
            plot((0:length(M(:,ll))-1)./Fs, M(:,ll),'color',[0,0,0]+1-(alpha^(k-ll))^2,'LineWidth',1.2);hold on;
        end
    else 
        %% 
        for ll = k-nb_frames_fading:k-1
            plot((0:length(M(:,ll))-1)./Fs, M(:,ll),'color',[0,0,0]+1-(alpha^(k-ll))^2,'LineWidth',1.2);hold on;
        end 
    end
    % Current display
    plot((0:length(M(:,k))-1)./Fs, M(:,k),'color',[0,0,0],'LineWidth',1.2);hold on %trace
    ylim([y_min_trace y_max_trace]); %Fix y-limits
    title(strcat('Exp: ',handles.experiment,'|Rec nb:',num2str(handles.recording_number),...
        '|Channel:',num2str(handles.channel),'|Cond:',num2str(handles.cond),'|Rep:',num2str(handles.repetition)))
    hold off;
    xl = xticklabels;
    xt = xticks;
    
    subplot(212);
    for i=1:size(frame_stim,1)
    plot(frame_time, frame_stim(i,:)); ylim([min_stim max_stim]); %Fix y-limits
    yticks([])
    hold on;  %stim
    end
    xticks([])
    %set(sp,'xtickl',[])
    %set(sp,'xticklabel',[])
    hold off;
    axis tight
    %xticks(xt), xticklabels(xl),
    xlabel('time[s]');
    lg = legend(handles.DAchan, 'Location','southeast');
    lg.FontSize = 8;
    


    % Save complete frame for video
    FV(k) = getframe(gcf);

end

% Parameterframe at the end of the movie

dim = [0.2 0.5 0.3 0.3];
str = {'Used movie parameters:',...
    strcat('BP filter: ','[',num2str(handles.filter.lowerthreshold),':',num2str(handles.filter.upperthreshold),']',' [Hz]'),...
    strcat('Peak TH: ',num2str(handles.peak.threshold),' [V]'),...
    strcat('Window size: ',num2str(handles.window_size.time),' [s]'),...
    strcat('DSR: ',num2str(display_speed_ratio),' [-]'),...
    strcat('Window step: ',num2str(handles.window_step_size.time),' [s]'),...
    strcat('Min displ time: ',num2str(handles.minimal_display_time.time),' [s]')};
annotation('textbox',dim,'Color','red','String',str,'BackgroundColor', 'k','FitBoxToText','on');
FV(k-1) = getframe(framesfig);
if strcmpi(framesfig.Visible,'off')
    framesfig.Visible = 'on';
    delete(framesfig);
end


video_frames = FV;
audio_frames_stim = FA_stim;
audio_frames_trace = FA_trace;
end

function close_showoff(src,callbackdata)
% Close request function 
% to hide figure and display a question dialog box 
% Close request function 
% to display a question dialog box 
   selection = questdlg('Close This Figure? It will interrump current set of frames.',...
      'Close frames figure',...
      'Yes','No','Yes'); 
   switch selection, 
      case 'Yes',
         delete(gcf)
      case 'No'
      return 
   end
end

