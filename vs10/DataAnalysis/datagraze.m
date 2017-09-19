function datagraze(Exp,varargin)
% datagraze - online analysis of recordings
%    datagraze('L15078b') opens a figure that shows online analysis 
%    for the recordings of the current dataset D of experiment L15078b,
%    i.e. L15078b_DS00NaN.EarlyDS.
%    The online analysis is performed every time the above file is changed.
%    The analyses performed depend on the information in D.
%
%    See also databrowse, savedataprogress.


% Set the matlab window title
addpath(genpath('C:\Early\vs10'),'-end');
rmpath(fullfile(matlabroot, 'toolbox', 'shared', 'statslib')); % removes dataset.m constructor in shared/statslib
clc;
disp(['--------------------------------------------------------------------' ...
    '-----------------------------------------------']);
disp('                                                   datagraze');
disp(['--------------------------------------------------------------------' ...
    '-----------------------------------------------']);
disp('');

try
    window_title = ['Datagraze - Experiment: ' Exp];
    jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    jDesktop.getMainFrame.setTitle(window_title);
end

if exist(Exp,'dir') % remote analysis
    Ddir = Exp;
else
    dummy_exp = experiment();
    Ddir = locate(dummy_exp, Exp); % local analysis
end
file = System.IO.FileSystemWatcher(Ddir); % look for changed files here
file.Filter = ['*_DS' zeropadint2str(nan,5) '.EarlyDS'];
file.EnableRaisingEvents = true;
file.NotifyFilter = System.IO.NotifyFilters.LastWrite;
addlistener(file,'Changed',@(Src,Ev,LR)local_doit(Src,Exp));
%h = msgbox({['Waiting for online analysis for experiment ' Exp '.'],...
%    'Use Ctrl-C to end.'},'Start');
while true
    pause(.5); % keep looking for changes
end

function local_doit(Src,Exp)
% The actual generation of the figure containing the data analyses.
%   Note: the "etime(Tlast,clock) < x" and "pause(x)" parameters were 
%   carefully chosen to perform online analysis as fast as possible,
%   but without running into concurrency issues.
persistent Tlast
persistent curTab
persistent i0 
persistent i1 
persistent figh
persistent analog_figh
persistent start_sample

%Src.EnableRaisingEvents = false; %disable filesystemwatcher

if ~isempty(Tlast) && (etime(clock, Tlast) < 2), return; end % ignore two or more related evetns
Tlast = clock;

if isempty(i0), i0 = 1; end
if i0 == 0, i0 = 1; end
if isempty(i1), i1 = inf; end
if isempty(figh), figh = -1; end
if isempty(analog_figh), analog_figh=-1; end
if isempty(start_sample), start_sample=1; end

has_analog_data = 'none';

try 
    pause(1);
    warning('off','MATLAB:cellRefFromNonCell');
%     warning ('off','all');
    D = read(dataset(), Exp, nan); 
%     warning ('on','all');
    warning('off','MATLAB:cellRefFromNonCell');
    if isfield(D.Stim,'DataGrazeActive')
        if D.Stim.DataGrazeActive == 1
            show_figures = 1;
        else
            show_figures = 0;
        end
    end
    
    if isfield(D.Data,'RX6_analog_1') && isfield(D.Data,'RX6_analog_2')
        has_analog_data = ['RX6_analog_1' ' & ' 'RX6_analog_2'];
        
    elseif isfield(D.Data,'RX6_analog_1') 
        has_analog_data = 'RX6_analog_1';
        
    elseif isfield(D.Data,'RX6_analog_2')
        has_analog_data = 'RX6_analog_2'; 
    end
    
    warning('on','MATLAB:cellRefFromNonCell');
catch
    D = dataset();
    disp('can''t read Dataset');
end;

if isvoid(D),  % in case no 00NaN DS was saved or it has been deleted
    disp(['Error while reading current dataset of experiment ', Exp,...
        ' for online analysis.']);
    return; 
end;

DVlist = cellify(D.Stim.Dataviewer); % make sure this is a cell
Pfile = cellify(D.Stim.DataviewerParamfile);

if show_figures == 1
    disp(['Current dataset of experiment ' Exp ' read.']);

    % Put everything together in a tabbed pane
    
    if figh == -1
        figh = figure;
    else
        figure(figh);
    end
    clf;
    placefig(figh,mfilename);

    warning('off','MATLAB:uitabgroup:OldVersion');
    warning('off','MATLAB:uitab:DeprecatedFunction');
    hTabGroup = uitabgroup('parent',figh, 'Units','normalized','Position',[0 0 1 1]); 
    % hTabGroup = uitabgroup('v0','parent',figh,'Units','normalized','Position',[0 0 1 1],...
    %     'SelectionChangeFcn', @(Src,Ev)local_selectionChangedFcn(Src,Ev,D,Pfile)); 
    drawnow;

    if isempty(curTab), temp = 1; else, temp = curTab; end; % store copy 

    for i=1:numel(DVlist)
        tab{i} = uitab(hTabGroup,'title',DVlist{i}); 
        eval([DVlist{i} '(D, tab{i}, Pfile{i});']);
        userdata_cell{i}=tab{i};
    end
    warning('on','MATLAB:uitab:DeprecatedFunction');
    warning('on','MATLAB:uitabgroup:OldVersion');
    % get the analog data
    if (strcmp(has_analog_data,'RX6_analog_1') || strcmp(has_analog_data,'RX6_analog_2') || ...
            strcmp(has_analog_data,'RX6_analog_1 & RX6_analog_2')) && strcmp(D.Stim.Analog_visualization,'active')
        
        % print the indices for debug purposes
%         fprintf('i0 = %d\t & i1 =%d\n',i0,i1);
        
        % Get the samples
        if strcmp(has_analog_data,'RX6_analog_1')
            [Y, dt, t0, Dtp] = samples(D.Data.RX6_analog_1, i0, i1);
            fs = Fsam(D.Data.RX6_analog_1);
        elseif strcmp(has_analog_data,'RX6_analog_2')
            [Y, dt, t0, Dtp] = samples(D.Data.RX6_analog_2, i0, i1);
            fs = Fsam(D.Data.RX6_analog_2);
        else
            [Y_tmp, dt, t0, Dtp] = samples(D.Data.RX6_analog_1, i0, i1);
            [Y_tmp2, dt, t0, Dtp] = samples(D.Data.RX6_analog_2, i0, i1);
            
            % in case of mismatch between the two lengths of datasets, pad
            % the longest with NaNs
            if length(Y_tmp) > length(Y_tmp2)
                length_diff = length(Y_tmp)-length(Y_tmp2);
                Y_tmp2(end+1:end+length_diff) = NaN;
            elseif length(Y_tmp) < length(Y_tmp2)
                length_diff = length(Y_tmp2)-length(Y_tmp);
                Y_tmp(end+1:end+length_diff) = NaN;
            end
            Y(1,:) = Y_tmp;
            Y(2,:) = Y_tmp2;
            fs = Fsam(D.Data.RX6_analog_1);
        end
%         % Do some moving average filtering
%         Y_orig = Y(5:end);
%         b = ones(5);
%         b = b/sum(b);
%         a =1;
%         Y = filter(b,a,Y);
%         Y=Y(5:end);
        
        % Update the indices
        i0 = length(Y); i1 = inf;
        

        
        
        % compute the time axis & plot the samples
        t = [0:length(Y)-1]*dt+t0;
        if analog_figh == -1
            analog_figh = figure;
        end
        figure(analog_figh);
        %set(analog_figh,'units','normalized','outerposition',[0 0 1 1]);
        clf;
        
        % Compute how many subplots are needed 1plot/2sec
        nchan = length(Y(:,1));
        colors ={'b','r'};
        for ichan=1:nchan
            subplot_interval = 2*fs; % samples
            Y_end = length(Y);
            num_plots = double(idivide(int32(Y_end),int32(subplot_interval),'ceil'));
            num_samples_last_plot = rem(Y_end,subplot_interval);
            if num_samples_last_plot == 0
                num_samples_last_plot = subplot_interval;
            end;
%             fprintf('num_plots=%d\n',num_plots);
            plot_index = 1+(ichan==2);
            for iplot=1:num_plots
%                 fprintf('iplot=%d\n',iplot);
                if iplot == num_plots
                    if num_plots == 1
                        plot(t,Y(ichan,:),colors{ichan});
                    else
                        % select the samples from Y
                        Y_sub = Y(ichan,int32(subplot_interval*(iplot-1)):end);

                        % Pad Y_sub with NaN's so the time axis of all subplots
                        % are as long 
                        num_nans = subplot_interval - length(Y_sub);
                        Y_sub(end+1:end+int32(num_nans)) = NaN;
                        % construct the time axis
                        time_offset = subplot_interval*(iplot-1)/fs*1000;
                        time_end = length(Y_sub);
                        t = ([0:time_end-1])*dt+t0+time_offset;

                        % create the subplot
                        subplot(num_plots,nchan,plot_index);
                        plot_index = plot_index+1+(nchan==2);
%                         fprintf('length of t=%d\t length of Y=%d\n',length(t),length(Y_sub));
                        plot(t,Y_sub,colors{ichan});
                    end
                else
                    % select the samples from Y
                    Y_sub = Y(ichan,int32(subplot_interval*(iplot-1)+1):int32(subplot_interval*iplot));

                    % construct the time axis
                    time_offset = subplot_interval*(iplot-1)/fs*1000;
                    time_end = length(Y_sub);
                    t = ([0:time_end-1])*dt+t0+time_offset;
                    % create the subplot
                    subplot(num_plots,nchan,plot_index);
                    plot_index = plot_index+1+(nchan==2);
%                     fprintf('length of t=%d\t length of Y=%d\n',length(t),length(Y_sub));
                    plot(t,Y_sub,colors{ichan});
                end
            end;
            xlabel('time in mili-seconds');
        end
        %plot(t,Y,'b',t,Y_orig,'r');
        
%         % Get the data splitted by icond per irep
%         fprintf('start_sample=%d   \t&& end_sample=%d\n',start_sample,start_sample+length(Y));
%         splitted_data = get_data_splits(D,D.Data.RX6_analog_1,start_sample, ...
%             start_sample + length(Y),Fsam(D.Data.RX6_analog_1),Y);
%         start_sample = start_sample + length(Y);

    end

    curTab = temp; % restore and use now the different tabs have been defined
    set(hTabGroup, 'SelectedIndex', curTab);
    set(hTabGroup, 'UserData', curTab);
    set(hTabGroup,'SelectionChangeCallback', @(Src,Ev)local_selectionChangedFcn(Src,Ev,D,Pfile,userdata_cell));

    %Src.EnableRaisingEvents = true; %enable filesystemwatcher
end

function local_selectionChangedFcn(Src,Ev,Dataset,Paramfile,UserDataCell)
figh = get(Src, 'parent');
itab = Ev.NewValue;
itab = find(UserDataCell == itab);
itab = itab(1);
disp(itab.Tag);
enableparamedit(Dataset, Paramfile{itab}, figh);
evalin('caller',['curTab = ', num2str(itab), ';']); % over different calls


function splitted_data = get_data_splits(ds,analog_data,start_sample, end_sample,Fsam,data)

splitted_data = [];
start_time = start_sample/Fsam*1000;
end_time = end_sample/Fsam*1000;

PRES = ds.Stim.Presentation;
onsets = PRES.PresOnset(2:end); % 2 as start index because 1 is baseline

[x,pres_ind] = NthFloor([start_time end_time],onsets);

icond_arr = ds.Stim.Presentation.iCond;
irep_arr = ds.Stim.Presentation.iRep;

ncond = PRES.Ncond;
nrep = PRES.Nrep;

for ip=1:PRES.Ncond*PRES.Nrep,       
    if ip == pres_ind(1)+1
        icond_start = icond_arr(ip+1);
        irep_start = irep_arr(ip+1);
    end
    
    if ip == pres_ind(2)
        icond_end = icond_arr(ip);
        irep_end = irep_arr(ip+1);
    end
    
end

% do some checking... Zzzz
if irep_start > irep_end
    error('irep_start > irep_end => wrong data indices');
end

if icond_start > ncond || icond_end > ncond
    error('icond_start > ncond of icond_end > ncond! icond_start: %d icond_end: %d', ...
        icond_start, incond_end);
end

if irep_start > nrep || irep_end > nrep
    error('irep_start > nrep of irep_end > nrep! irep_start: %d irep_end: %d', ...
        irep_start, inrep_end);
end
% iterate over each rep

% iterate in irep_start from icond+1 (if icond+1<ncond) to ncond
% iterate in irep_end from icond=1 to icond=icond_end
% iterate in other irep from icond=1 to icond=ncond

for irep=irep_start:irep_end
   
    if irep == irep_start

        % first icond is probably not completed so start from icond+1 if
        % possible
        if icond_start+1 <= ncond
            icond_start=icond_start+1;
        end
        
        % Check if nrep_start == nrep_end if so icond should be limited to
        % icond_end
        if irep_start == irep_end
           ncond = icond_end; 
        end
        
        % loop over all the conditions
        for icond=icond_start:ncond
            % debug
%             fprintf('icond=%d   \tirep=%d\n',icond,irep);
            Pres_onset = round(onsets(icond*irep)*Fsam/1000);
            Pres_end = round(onsets(icond*irep+1)*Fsam/1000);
%             fprintf('length of data buffer=%d   \tlower bound=%d   \t upper bound=%d\n', ...
%                 length(data),int32(Pres_onset), int32(Pres_end));
            splitted_data{irep,icond} = data(int32(Pres_onset):int32(Pres_end)); 
        end
        
    elseif irep == irep_end
        % last icond is probably not completed so end at icond-1 if
        % possible
        if icond_end-1 > 0
            icond_end=icond_end-1;
        end
        
        % loop over all the conditions
        for icond=1:icond_end
            % debug
%             fprintf('icond=%d   \tirep=%d\n',icond,irep);
            Pres_onset = round(onsets(icond*irep)*Fsam);
            Pres_end = round(onsets(icond*irep+1)*Fsam);
%             fprintf('length of data buffer=%d   \tlower bound=%d   \t upper bound=%d\n', ...
%                 length(data),int32(Pres_onset), int32(Pres_end));
            splitted_data{irep,icond} = data(Pres_onset:Pres_end); 
        end
    else
        % Normally all icond's for this irep are available for this irep
        for icond=1:ncond
            % debug
%             fprintf('icond=%d   \tirep=%d\n',icond,irep);
            Pres_onset = round(onsets(icond*irep)*Fsam);
            Pres_end = round(onsets(icond*irep+1)*Fsam);
%             fprintf('length of data buffer=%d   \tlower bound=%d   \t upper bound=%d\n', ...
%                 length(data),int32(Pres_onset), int32(Pres_end));
            splitted_data{irep,icond} = data(Pres_onset:Pres_end); 
        end
    end
    
end


