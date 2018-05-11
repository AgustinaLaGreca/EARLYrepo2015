function [] = window_size_and_plot(hObject,handles)
% [] = window_size_and_plot(hObject,handles)
%   window_size_and_plot plots the result of the selected windowsize (handles.window_size)


%% Window size

% Window size should be larger than 1 sample (time) = 1/Fs
if handles.window_size.samples < 1
    msgbox(strcat('The window size should be larger than one time step (=',num2str(1/handles.Fs),'s), please change your input.'));
    return
end
% Window size should be smaller than the total signal (without margin) =
% length(trace)*Fs - window size
if handles.window_size.samples > handles.nb_of_samples/2
    msgbox(strcat('The window size should be smaller than half the total signal (=',num2str(handles.nb_of_samples/2/handles.Fs),'s), please change your input.'));                                                             
    return
end

%% -> Selection of the window size

    window_sizes = handles.window_size.samples; %from s to nb of samples
    preplot_peak_locs = handles.peak.positions(handles.peak.positions > window_sizes/2);
    preplot_peak_locs = preplot_peak_locs(preplot_peak_locs < handles.nb_of_samples-window_sizes/2);

    nb_of_plots = min([7,length(preplot_peak_locs)]);
    if isempty(handles.peak.positions)
            msgbox(strcat('No peaks are selected (by the current threshold).'));
    elseif nb_of_plots == 0 
            msgbox(strcat('Window to large to plot a window, please reduce windowsize.'));

    else
        figure();
%         ylim = max(handles.trace);
%         plot(handles.t,handles.trace); hold on
        
%         C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};
%         for ii = 1:nb_of_plots
%             frame_ind_lower = (preplot_peak_locs(ii)-floor((window_sizes)/2))/handles.Fs;
%             frame_ind_higher = (preplot_peak_locs(ii)+floor((window_sizes-1)/2))/handles.Fs;
%             line([frame_ind_lower frame_ind_lower], ylim,'color',C{ii});
%             line([frame_ind_higher frame_ind_higher], ylim,'color',C{ii});
%             axis tight
%         end
%         hold off
            for i = 1:nb_of_plots
            frame_ind = ((preplot_peak_locs(i)-floor((window_sizes)/2)):1:...
                         (preplot_peak_locs(i)+floor((window_sizes-1)/2)));
            subplot(nb_of_plots,1,i)
            % from 2016b
%             plot((1:window_sizes)./handles.Fs,handles.trace(frame_ind),...
%                         '-*','MarkerIndices',floor((window_sizes)/2)+1,'MarkerEdgeColor','r')
            plot((1:window_sizes)./handles.Fs,handles.trace(frame_ind)) %before 2016b
            axis tight
            xlabel('time [s]')

            end
    end

end
