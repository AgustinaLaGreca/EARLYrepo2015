function [peak_locs, stim_int, t_int] = peak_finder_and_plot(hObject,handles,plot_demand)
% [peak_locs] = peak_finder_and_plot(hObject,handles,plot_demand)
%   peak_finder_and_plot finds the peak values in hadles.trace above a
%   certain threshold (handles.peak.threshold) and returns the location of all the found
%   peaks in an array (peak_locs). If plot_demand=1, the result of
%   the peak search action is displayed on the gui. If plot_demand=0, nothing is
%   displayed.

% Peaks are found in the trigger channel selected. The structure of peaks
% has the following fields:
% handles.peak.threshold value for triggering
% handles.peak.channel : int specifyng the channel to trigger:
%       channel = 1 is triggering on the trace, channel = {2,3} is
%       triggering on handles.stim(channel-1,:)
% handles.peak.positions array 

%% Peak detection asserts

% Threshold should be positive
if handles.peak.threshold < 0
    msgbox('Please insert a positive value for the threshold, please change your input.');
    peak_locs = [];
    return
end

%% Triggering channel
if handles.peak.channel > 1
    trigChannel = handles.stim(handles.peak.channel-1,:);
    peak_locs = findtrigger(trigChannel, handles.peak.threshold, 'positive');
    
    % Interpolate channel, time and other stimulus if exists
    t_int_vals = interpolatechannel(trigChannel, handles.t, peak_locs, handles.peak.threshold);
    t_int = handles.t;
    t_int(peak_locs) = t_int_vals;
    
    stim_int = handles.stim;
    if size(handles.stim,1) > 1
        otherChannel = 1:size(handles.stim,1) ~= (handles.peak.channel-1);
        otherChannel = find(otherChannel);
        for i=otherChannel
            ochannel_int_vals = interpolatechannel(handles.t, handles.stim(i,:), peak_locs, t_int_vals);
            stim_int(i,peak_locs) = ochannel_int_vals;
        end
    end
    stim_int(handles.peak.channel-1,peak_locs) = handles.peak.threshold;
    
else
    %% Triggering on trace: finding peaks
    trigChannel = handles.trace;
    % Preview of the whole signal and selected peaks
    % peakfinder
    warning('off','signal:findpeaks:largeMinPeakHeight') %turn the warning off in case no peaks are found
                                                       % (threshold too
                                                       % high)
    [~,peak_locs] = findpeaks(trigChannel,'MinPeakHeight',handles.peak.threshold);
    stim_int = [];
    t_int = [];
    
end
%[~,peak_locs] = findpeaks(trigChannel,'MinPeakHeight',handles.peak.threshold);

%% Plot the data
if plot_demand
    
    % Y limits trace
    if handles.y_limits_trace.auto_on
        y_min_trace = handles.y_limits_trace.auto.min;
        y_max_trace = handles.y_limits_trace.auto.max;
    else
        y_min_trace = handles.y_limits_trace.manual.min;
        y_max_trace = handles.y_limits_trace.manual.max;
    end
    
%   plot(handles.t,handles.trace,'-*','MarkerIndices',peak_locs,'MarkerEdgeColor','r');hold on %from2016b

    % alternative for before 2016b
    MarkerIndices = peak_locs;  
    plot(handles.t,handles.trace, 'b-'); hold on;               %plot everything with appropriate line color and no marker
    plot(handles.t(MarkerIndices), handles.trace(MarkerIndices), 'r*');  %plot selectively with appropriate color and marker but no line
    axis tight
    
    ylim([y_min_trace y_max_trace]); %Fix y-limits
    try
    hline = refline(0,handles.peak.threshold); 
    hline.Color = 'g';
    catch ME
        disp(['Threshold line not displayed. Error: ', ME.message]);
    end
    hold off;
    xlabel('time [s]')
end

end

function indexes = findtrigger(yvalues, threshold, varargin)
%% FINDTRIGGER
% Return indexes of values vector that would work as trigger points. Those
% points are the points equal to or the colsest inmediately after that
% threshold.
% yvalues are the signal y values
% threshold value for trigger
% Varargin is the POLARITY, 'positive' or negative'. 'positive' as default

if nargin == 2
    polarity = 'positive';
else
    polarity = varargin{1};
end
higherInd = [0 diff(yvalues > threshold)];
if strcmpi(polarity, 'positive')
    indexes = find(higherInd > 0);
else
    indexes = find(circshift((higherInd < 0),numel(higherInd)-1));
end

end

function xvals = interpolatechannel(yvalues, tvalues, trig_indexes, value)
% Interpolate a tvalue using yvalues, tvalues around trig_indexes.
% The value interpolated is value. It can be a single value or an array of
% the same size as trig_indexes

intvalues = 2;
xvals = zeros(size(trig_indexes));
if size(value) == 1
    value = value*ones(size(trig_indexes));
end
for j=1:numel(trig_indexes)
    i = trig_indexes(j);
    ind = i-intvalues:i+intvalues;
    xinter = interp1(yvalues(ind),tvalues(ind),value(j),'pchip');
    xvals(j) = xinter;
end
end
