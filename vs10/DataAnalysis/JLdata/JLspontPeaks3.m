function [ST, Y] = JLspontPeaks3(Uidx, ThrPeak, MinDT, doPlot);
% JLspontPeaks3 - counting and analyzing peaks of spontaneous activity
%   Y=JLspontPeaks2(Uidx, DTsmooth, MinDist_ms)

Fhighpass = 100; % Hz highpass filtering
t_stimOnset = 500; % ms stim onset
DTsmooth = 0.08; % ms

% default input args
[ThrPeak, MinDT, doPlot] = arginDefaults('ThrPeak/MinDT/doPlot', [], [], 0);

CD = 'D:\processed_data\JL\spontEPSP'; % cache dir
CFN = fullfile(CD, [mfilename num2str(Uidx)]);
[Y, CFN, CP] = getcache(CFN, {Fhighpass, MinDT,  ThrPeak, DTsmooth});

if isempty(Y) || doPlot, % compute it
    % AP criterium
    [S, rec] = JLgetRec(Uidx);
    APcrit = [S.APthrSlope S.APwindow_start S.APwindow_end S.CutoutThrSlope];
    if numel(APcrit)<4, APcrit(4) = APcrit(1); end  % default trunc Thr is AP detection thr

    % Default ThrPeak, DTsmooth vay with cell. Get from local "database".
    [ThrPeak, MinDT] = local_param(S.icell_run, ThrPeak, MinDT);

    %read data, truncate spikes, highpass filter & smoothen
    dt = S.dt;
    rec = rec-mean(rec(:));
    Nsam = round(t_stimOnset/dt);
    [rec, tAP] = APtruncate2(rec(1:Nsam+100), APcrit(4), dt, APcrit(2:3));
    [B,A] = butter(5, 2e-3*Fhighpass*dt ,'high');
    rec = filter(B,A,rec);
    rec = smoothen(rec, DTsmooth, dt);
    rec = rec(1:Nsam); % stop @ stim onset

    % find peaks that are at least MinDT ms apart
    NsamDist = ceil(MinDT/dt); % min distance to next peak
    NsamDip = ceil(1.5*MinDT/dt); % range over which to evaluate size
    [RMX, ipeak] = runmax(rec, NsamDist); % detect peaks
    RMN = -runmax(-rec, NsamDip); % % effective size of peaks
    PeakSize = RMX(ipeak)-RMN(ipeak);
    qbigpeak = PeakSize>ThrPeak;
    Tpeak = dt*(ipeak-1);
    if doPlot,
        set(figure,'units', 'normalized', 'position', [0.082 0.482 0.901 0.417]);
        dplot(dt, rec, 'k');
        xdplot(dt, RMX, 'b:');
        xdplot(dt, RMN, 'b:');
        fenceplot(Tpeak, ylim, 'g');
        fenceplot(Tpeak(qbigpeak), ylim, 'r');
        TracePlotInterface(gcf);
        if doPlot==-1, return; end;
        waitfor(gcf);
    end

    % restrict to big peaks
    ipeak = ipeak(qbigpeak);
    Tpeak = dt*(ipeak-1);
    Vpeak = rec(ipeak);
    PeakSize = PeakSize(qbigpeak);
    Npeak = numel(ipeak);
    Nsam = numel(rec);

    % half widths
    HalfHeight = Vpeak-0.5*PeakSize;
    HalfWidth = nan(1,Npeak);
    for ip=1:Npeak,
        hh = HalfHeight(ip);
        isam = ipeak(ip);
        i0 = max(1,isam-NsamDip);
        i1 = min(Nsam,isam+NsamDip);
        ileft = find(rec(i0:isam)>hh, 1, 'first')+i0;
        irite = find(rec(isam:i1)<hh, 1, 'first')+isam;
        if ~isempty(ileft) && ~isempty(irite),
            HalfWidth(ip) = dt*(irite-ileft);
        elseif ~isempty(ileft),
            HalfWidth(ip) = 2*dt*(isam-ileft);
        elseif ~isempty(irite),
            HalfWidth(ip) = 2*dt*(irite-isam);
        end
    end
    % reject excessively broad peaks
    qok = HalfWidth<=2*MinDT;
    %dsize(qok, PeakSize, HalfWidth, Tpeak);
    PeakSize = PeakSize(qok); HalfWidth = HalfWidth(qok); Tpeak = Tpeak(qok);
    Npeak = numel(Tpeak);
    [SpontEPSP_DTsmooth, SpontEPSP_Fhighpass, SpontEPSP_MinDT] = deal(DTsmooth, Fhighpass, MinDT);
    UniqueRecordingIndex = Uidx;
    NspontEPSP = Npeak;
    SpontEPSPsize = PeakSize(:);
    EPSPhalfWidth = HalfWidth(:);
    TspontEPSP = Tpeak(:);
    icell_run = S.icell_run;
    Y = CollectInStruct(UniqueRecordingIndex, icell_run, '-SpontEPSPs', ...
        ThrPeak, SpontEPSP_DTsmooth, SpontEPSP_Fhighpass, SpontEPSP_MinDT, ...
        SpontEPSPsize, EPSPhalfWidth, NspontEPSP, TspontEPSP);
    putcache(CFN, 1, CP, Y);
end

UniqueRecordingIndex = Uidx;
icell_run = Y.icell_run;
NspontEPSP = Y.NspontEPSP;
meanSpontEPSPsize = mean(Y.SpontEPSPsize);
stdSpontEPSPsize = std(Y.SpontEPSPsize);
meanEPSPhalfWidth = mean(Y.EPSPhalfWidth);
stdEPSPhalfWidth = std(Y.EPSPhalfWidth);

ST = CollectInStruct(UniqueRecordingIndex, icell_run, '-SpontEPSPs', ...
    NspontEPSP, meanSpontEPSPsize, stdSpontEPSPsize, meanEPSPhalfWidth, stdEPSPhalfWidth);
DBFN = fullfile(CD, 'JLspontPeaks.dbase');
if ~exist(DBFN, 'file'),
    init_dbase(DBFN, 'UniqueRecordingIndex');
end
Add2dbase(DBFN, ST);

%==============================================
function [ThrPeak, MinDT] = local_param(icell_run, ThrPeak, MinDT);
if isempty(MinDT), MinDT = 0.3; end % % ms min distance btw EPSPs
if isempty(ThrPeak), % get from table
    ThrPeak = 0.2;
    switch icell_run,
        case 1,  ThrPeak = 0.15;
        case 4,  ThrPeak = 0.15; MinDT = 0.35;
        case 30,  ThrPeak = 0.2;
        otherwise,
            warning(['No default ThrPeak known for this cell: icell_run = ' num2str(icell_run)]);
    end
end

