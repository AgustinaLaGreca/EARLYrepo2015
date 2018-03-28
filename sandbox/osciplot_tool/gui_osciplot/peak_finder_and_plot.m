function [peak_locs] = peak_finder_and_plot(hObject,handles,plot_demand)
% [peak_locs] = peak_finder_and_plot(hObject,handles,plot_demand)
%   peak_finder_and_plot finds the peak values in hadles.trace above a
%   certain threshold (handles.peak.threshold) and returs the location all the found
%   paeks in an array (peak_locs). If plot_demand=1, the result of
%   the peak search action is displayed on the gui. If plot_demand=0, nothing is
%   displayed.

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
if plot_demand
    
%   plot(handles.t,handles.trace,'-*','MarkerIndices',peak_locs,'MarkerEdgeColor','r');hold on %from2016b

    % alternative for before 2016b
    MarkerIndices = peak_locs;  
    plot(handles.t,handles.trace, 'b-');   hold on               %plot everything with appropriate line color and no marker
    plot(handles.t(MarkerIndices), handles.trace(MarkerIndices), 'r*');  %plot selectively with appropriate color and marker but no line
    axis tight
    hline = refline(0,handles.peak.threshold); 
    hline.Color = 'g';hold off
    xlabel('time [t]')

end
end


