function [peak_locs] = peak_finder_and_plot(hObject,handles);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Peak detection asserts

% Threshold should be positive
if handles.peak.threshold < 0
    msgbox('Please insert a positive value for the threshold, please change your input.');
    
    peak_locs = [];
    return
end

%% Preview of the whole signal and selected peaks

    %peakfinder
    warning('off','signal:findpeaks:largeMinPeakHeight') %turn the warning off in case no peaks are found
                                                       % (threshold to
                                                       % high)
    [~,peak_locs] = findpeaks(handles.trace,'MinPeakHeight',handles.peak.threshold);

    % Plot the data
    plot(handles.t,handles.trace,'-*','MarkerIndices',peak_locs,'MarkerEdgeColor','r');hold on
    hline = refline(0,handles.peak.threshold); 
    hline.Color = 'g';hold off
%     axis tight
    xlabel('time [t]')


end


