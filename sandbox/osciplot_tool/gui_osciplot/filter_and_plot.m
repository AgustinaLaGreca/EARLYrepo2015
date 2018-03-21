function [trace_filtered] = filter_and_plot_filter(hObject,handles)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


%% filter asserts

% LL and UL should be positive
if handles.filter.lowerthreshold <= 0 || handles.filter.upperthreshold <= 0
    msgbox('Please insert strictly positive values for the cuttof frequencies of the filter, please change your input.');
    
    trace_filtered = handles.trace;
    return
end

% LL should be smaller than UL
if handles.filter.lowerthreshold >= handles.filter.upperthreshold 
    msgbox('The lower CUF should be smaller than the upper CUF, please change your input.');
    
    trace_filtered = handles.trace;
    return
end


% UL should be smaller than the Nyquist frequency (=0.5*Fs)
if handles.filter.upperthreshold >= handles.Fs/2
    msgbox(strcat('The upper CUF should be smaller than the Nyquist freq (=',num2str(floor(handles.Fs/2)),'Hz), please change your input.'));
    
    trace_filtered = handles.trace;
    return
end

%% filter itself

    % Bandpass Butterworth fiter
    cof =  [handles.filter.lowerthreshold,handles.filter.upperthreshold]; %cut off freq
    [b,a] = butter(2,cof./handles.Fs);
    trace_filtered = filtfilt(b,a,handles.trace_original);
    
    % result of filtering
    plot(handles.t,handles.trace_original);hold on
    plot(handles.t,trace_filtered);hold off
    axis tight
    xlabel('time [t]')
    legend('original trace','filtered trace')
    
end


