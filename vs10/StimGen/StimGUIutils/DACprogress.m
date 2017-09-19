function Y=DACprogress(keyword, X, figh);
%  DACprogress - report progress of D/A conversion
%   Preparation calls:
%     DACprogress('feedback_fcn', @foo) stores the feedback function which
%     provides info on stimulus conditions. This call is typically issued
%     by sortConditions, but the assigment of the feedback function may be
%     overruled by stimulus generators makestimXXX.
%
%     DACprogress('condition_array', iCond) stores array iCond. iCond(k) is
%     the condition index belonging to the kth presentation. See
%     sortConditions.
%
%     DACprogress('nsam_info', Nsam, figh) stores the number of samples of
%       each stimulus condition. Nsam is an array. Nsam(k) is the # samples
%       making up the kth stimulus condition. The index k refers to sored
%       waveforms, not to the order of play. Figh is handle of the GUI to 
%       which stimulus info will be sent. This type of call is issued by
%       Waveform/play.
%
%   Action calls:
%    DACprogress('report') or DACprogress calls seqplayStatus to know which
%    buffer is played, uses the buffer info to convert this to the stimulus
%    index, and uses the feedback function to display info on this stimulus
%    condition. The info is displayed to the GUI with handle figh, or to
%    the Matlab command window if figh is empty. This type of call is
%    not issued directly, but only ondirectly via the 'start_timer' call.
%
%    DACprogress('start_timer', DT)  starts a Matlab timer which reports
%    the D/A progress every DT ms by calling DACprogress('report').
%    The default interval DT is 100 ms. When the D/A has stopped, this is 
%    reported, and the timer is stopped.
%
%    DACprogress('stop') aborts the current, timer-based, reporting and
%    gives a message "D/A interrupted".
%    DACprogress('stop', Str) replaces the standard warning by Str.
%    DACprogress('stop', Str, MessMode) also controls way in which the 
%    Str message is displayed. For valid disply Mode, see GUImessage.
%
%   Clear:
%     DACprogress('clear') clears all previously stored info.
%
%    See also SortConditions, Waveform/play, GUImessage.

persistent DACinfo
%'DACprogress TMP disabled', return
%'DACprogress TMP enabled'
if nargin<1, keyword = 'report';
end
switch lower(keyword),
    case 'init', % initialize by storing to two functions that enable feedback
        DACinfo = [];
        DACinfo.iSam2iPres = X.iSam2iPres; % DAC sample count --> stim presentation count
        DACinfo.DisplayFcn = X.DisplayFcn; % stim presentation count --> feedback string
        DACinfo.figh = figh; 
    case 'clear',
        DACinfo = [];
    case 'debug',
        Y = DACinfo;
    case 'report',
        Y = [];
        return;
        S = seqplaystatus;
        if ~S.Active, Str = 'No D/A.';
        else,
            RemTime = (S.NsamTot-S.isam_abs)*S.dt; % ms remaining DA time
            if isempty(DACinfo),
                Str{1,1} = ['Playing from Wavebuffer ' num2str(S.WaveIndex) '.'];
            else,
                ipres = DACinfo.iSam2iPres(S.isam_abs);
                Str{1,1} = DACinfo.DisplayFcn(ipres);
            end
            Str{2,1} = ['Remaining D/A time: ' DurString(RemTime)];
            Y = RemTime;
        end
        local_disp(Str, DACinfo.figh);
    case 'start_timer',
        if nargin<2, DT=nan; else, DT=X; end 
        if nargin<3, StopFcn=''; else, StopFcn=figh; end % no stopfcn for timer
        if isnan(DT), DT = 100; end % 100 ms default refresh interval
        DACinfo.Timer = timer('Name', 'DAcprogressReporter', 'TimerFcn','DACprogress exec_timer', ...
            'Period', DT/1e3, 'ExecutionMode', 'fixedSpacing', 'StopFcn', StopFcn);
        start(DACinfo.Timer);
    case 'exec_timer',
        TT=DACprogress('report');
        if isempty(TT),
            stop(DACinfo.Timer);
            delete(DACinfo.Timer);
            local_disp('D/A finished', DACinfo.figh);
        end
    case 'stop', % DACprogress('stop', Str, MessMode)
        if nargin<2, X = 'D/A interrupted'; end
        if nargin<3, figh = 'warning'; end
        Mess = X; % abuse of 2nd arg
        MessMode = figh; % abuse of 3rd arg
        if isvalid(DACinfo.Timer),
            stop(DACinfo.Timer);
            delete(DACinfo.Timer);
        end
        local_disp(Mess, DACinfo.figh, MessMode);
    otherwise,
        error(['Invalid keyword ''' keyword '''.']);
end

%==============================================
function local_disp(Str, figh, MessMode);
if nargin<3, MessMode='shy'; end % don't bring GUI to the fore all the time - we're..
if isempty(figh), disp(Str); % ... probably running from the background.
else, GUImessage(figh, Str, MessMode); 
end 





