function Y = JLspikeThr(Uidx, SlopeThr, figh);
% JLspikeThr - interactive assessment of spike thresholding
%    JLspikeThr(Uidx, SlopeThr)
%    SlopeThr in V/s is threshold slope of dV/dt.
%    If SlopeThr<0 (the usual case), APs are identified based on
%    dV/dt<=SlopeThr.
%    If SlopeThr>0, APs are identified based on dV/dt>=SlopeThr.
%    The default SlopeThr value is the value previously determined, which
%    is retrieved by the APthrSlope field of JLdatast(Uidx).
%    To use this defaulyt value, either do not specify the input arg
%    SlopeThr or use SlopeThr=nan.
%
%    Y = JLspikeThr(Uidx, SlopeThr, figh) plots to figure with handle figh. 
%    If figh==0, no plotting is done. The output arg Y is a struct with 
%    fields
%      UniqueRecordingIndex: 204148005
%                  SlopeThr: AP threshold  in V/s (equal to input SlopeThr)
%                  doInvert: true if SlopeThr is negative.
%                    PolFac: 1 or -1 for doInvert = false, true.
%                    minSlp: large collection of events used to compile
%                            histogram of slopes.
%                    Tevent: times (ms from start of recording) of a
%                            collection of largest events that spans ~equal
%                            amount of supra thresholds and subthreshold
%                            events, or a handful of subthreshold ones if
%                            the recording does not contain any
%                            suprathreshold events. Events are timed based
%                            on their peak just preceding the steepest part
%                            of the slope. Events are sorted from large to
%                            small steepst slopes.
%            EventDownSlope: Steepest slopes of the events occurring at 
%                            Tevent. Desending array due to sorting.
%                  t_offset: latencies between steepest part and peak of
%                            each of the events.
%                     Tsnip: singe time (ms) array for plotting event
%                            snippets.
%                   APsnips: Event snippets. APsnips(:,I) is the waveform
%                            of the event occurring at Tevent(I).
%                       qAP: logical array indicating which events exceed
%                            the threshold SlopeThr.
%                    Tdsnip: Time array (ms after rec start) for time
%                            derivative of waveforms.
%                  dAPsnips: snippets of time derivatives of events.
%
%    JLspikeThr(Uidx, nan, 'rmcache') removes the cache.
%
%
%    See also JLspikes, JLspikeStats, JLspikeAlign.

[SlopeThr, figh] = arginDefaults('SlopeThr/figh', nan, []);

Y = local_doit(Uidx, SlopeThr, figh);
if ~isequal(0,figh) && ~ischar(figh),
    local_plot(Y, figh);
end

function Y = local_doit(Uidx, SlopeThr, figh);
S = JLgetRec(Uidx); % this call is fast due to caching
if isnan(SlopeThr), SlopeThr = S.APthrSlope; end
doInvert = SlopeThr<0;
SlopeThr = abs(SlopeThr);

icell = getfield(JLdbase(Uidx), 'icell_run');
CFN = fullfile(processed_datadir, 'JL', 'JLspikethr', ['cell_' num2str(icell)]);
[Y, CFN, CP] = getcache(CFN, {Uidx doInvert});
if isequal(pi, Y), Y = []; end
if isequal('rmcache', figh),
    putcache(CFN, 300, CP, pi); % remove this cache
    return;
end
if isempty(Y), % isolate events using one third of the current default threshold 
    [S, rawrec] = JLgetRec(Uidx); % rawrec has no APs or trend removed. It's only 50-Hz filtered.
    drr = diff(smoothen(rawrec, 0.05, -S.dt)/S.dt); % time derivative
    [drr(1) drr(end)] = deal(drr(2), drr(end-1)); % correct boundary discontinuities
    if doInvert, PolFac = -1; else, PolFac=1; end
    [Tneg, minSlp, isort, RM] = peakfinder(S.dt, PolFac*drr, -0.8, PolFac*S.APthrSlope/3); % minus sign rejects closely spaced events
    [Tneg, minSlp] = deal(Tneg(isort), minSlp(isort)); % sort from big to small
    irm = find(Tneg<=S.dt);
    Tneg(irm)=[]; minSlp(irm)=[]; 
    %plot(minSlp);
    NAP = sum(minSlp>=SlopeThr); % # accepted APs
    %numel(Tneg), 2*NAP
    Ntot = min(numel(Tneg), max(10, 2*NAP)); % max # events to be shown below
    [Tevent, EventDownSlope] = deal(Tneg(1:Ntot), minSlp(1:Ntot));
    % isolate AP candidate, aligning their peaks
    [dum, Tsnip, dAPsnips] = cutFromWaveform(S.dt, drr, Tevent, [-2 0]);
    t_offset = 0; % avoid crash when Tevent is empty (no APs)
    for iAP=1:numel(Tevent),
        ipeak = find(dAPsnips(:,iAP)>0,1,'last'); % latest sample index at which potential was rising
        if isempty(ipeak), ipeak = 1; end % default for "peakless" segments; see a few lines down
        t_offset(iAP,1) = Tsnip(ipeak); % delay between AP peak and instant steepest downward slope
    end
    Tpeak = Tevent + t_offset;
    irm = find(Tpeak+S.APwindow_end<=0); % remove before range range segments, potentially resulting from peakless cases
    Tpeak(irm) = []; Tevent(irm) = []; t_offset(irm) = []; EventDownSlope(irm) = [];
    CutWin = [S.APwindow_start S.APwindow_end];
    [Yrem, Tsnip, APsnips] = cutFromWaveform(S.dt, rawrec, Tpeak, CutWin);
    [dum, Tdsnip, dAPsnips] = cutFromWaveform(S.dt, drr, Tpeak, CutWin);
    UniqueRecordingIndex = Uidx;
    qAP = nan;
    Y = CollectInStruct(UniqueRecordingIndex, '-AP_thresholding', SlopeThr, doInvert, PolFac, ...
        minSlp, Tevent, EventDownSlope, t_offset, Tsnip, APsnips, qAP, Tdsnip, dAPsnips);
    putcache(CFN, 300, CP, Y); % no cell has more than 290 recs
end
Y.SlopeThr = SlopeThr;
Y.qAP = Y.EventDownSlope>=SlopeThr; % accepted APs
%sum(qAP)
%dsize(Tpeak, Tsnip, APsnips ,qAP)
%dplot(S.dt, rawrec); traceplotInterface(gcf); fenceplot(Tpeak(qAP), ylim ,'g'); fenceplot(Tpeak(~qAP), ylim ,'r'); 


function figh = local_plot(Y, figh);
if isempty(figh), 
    figh = figure;
    set(figh,'units', 'normalized', 'position', [0.0727 0.333 0.825 0.589]);
end
Uidx = Y.UniqueRecordingIndex;
% ==histogram of slopes, working exclusively from raw slope stats in variable minSlp
logSlope = log(Y.minSlp(Y.minSlp>0));
Nbin = min(50, numel(logSlope)/10);
BC = linspace(min(logSlope), max(logSlope), Nbin);
subplot(2,4,1:2);
try,
    hist(logSlope, BC); NN = hist(logSlope, BC); % 1st call: force plot; 2nd: get hist values.
    xlim([min(BC) max([max(BC), 1.05*log(Y.SlopeThr)])]);
    Xval = [1 2 3 4 5 7 8 10 12 15 20 25 30 50 100];
    set(gca,'xtick', log(Xval), 'xticklabel', Words2cell(num2str(Xval)), 'fontsize',11);
    fenceplot(log(Y.SlopeThr), ylim, 'color', [0 0.75 0], 'linewidth',2);
    Peak2 = max(NN(BC>=log(Y.SlopeThr)));
    Ymax = min(max(NN), 10*Peak2);
    if ~isempty(Ymax), ylim([0 Ymax]); end
end
if Y.doInvert, xlabel('Abs downward slope (V/s)', 'fontsize', 12);
else, xlabel('Upward slope (V/s)', 'fontsize', 12);
end
ylabel('# occurrences');
% ==spike and nonspike waveforms
subplot(2,4,3:4);
% if isempty(Y.APsnips),
%     text(0.5, 0.5, 'No APs', 'fontsize', 14, 'color', 'r');
%     return;
% end
NplotMax = 20;
% plot at most NplotMax supra- and subthreshold events
[iabove, Nabove] = deal(find(Y.qAP), sum(Y.qAP));
iabove = iabove(max(1,Nabove-NplotMax+1):Nabove);
[ibelow, Nbelow] = deal(find(~Y.qAP), sum(~Y.qAP));
ibelow = ibelow(1:min(Nbelow,NplotMax));
if ~isempty(ibelow), plot(Y.Tsnip, Y.APsnips(:,ibelow), 'b'); end
if ~isempty(iabove), xplot(Y.Tsnip, Y.APsnips(:,iabove), 'r'); end
% avoid complete asking of non APs by APs: replot half of the nonAPs on top of the APS
if ~isempty(ibelow), xplot(Y.Tsnip, Y.APsnips(:,ibelow(1:2:end)), 'b'); end
ylabel('V (mV)', 'fontsize', 12)
text(0,1.1, JLtitle(Uidx), 'horizontalalign', 'center', 'units' ,'normalized', 'fontweight', 'bold', 'fontsize', 14);
subplot(2,4,7:8);
if ~isempty(ibelow), plot(Y.Tdsnip, Y.dAPsnips(:,ibelow), 'b'); end
if ~isempty(iabove), xplot(Y.Tdsnip, Y.dAPsnips(:,iabove), 'r'); end
% again, replot half of the nonAPs on top of the APS
if ~isempty(ibelow), xplot(Y.Tdsnip, Y.dAPsnips(:,ibelow(1:2:end)), 'b'); end
xlabel('Time (ms)', 'fontsize', 12)
ylabel('dV/dt (V/s)', 'fontsize', 12)
% update stuff
setGUIdata(figh, 'Uidx', Uidx);
uicontrol('style', 'text', 'position', [60 120 80 25], 'string', 'Threshold: ', 'fontsize', 12)
hEdit = uicontrol('style', 'edit', 'backgroundcolor', 'w' ,'position', [150 120 100 25], ...
    'fontsize', 12, 'horizontalalign', 'left', 'string', num2str(Y.PolFac*Y.SlopeThr));
uicontrol('style', 'text', 'position', [230 120 40 25], 'string', 'mV ', 'fontsize', 12)
uicontrol('style', 'pushbutton', 'position', [280 120 80 25], 'string', 'Replot', 'fontsize', 13 ,...
    'callback', @(Src,Ev)local_replot(hEdit, parentfigh(Src)));
uicontrol('style', 'pushbutton', 'position', [400 120 80 25], 'string', 'Empty cache', 'fontsize', 13 ,...
    'callback', @(Src,Ev)local_rmcache(hEdit, parentfigh(Src)));
uicontrol('style', 'pushbutton', 'position', [200 80 150 30], 'string', 'Update', 'fontsize', 13 ,...
    'BackgroundColor', [0.7 0.2 0.2], 'callback', @(Src,Ev)local_update(hEdit));
uicontrol('style', 'pushbutton', 'position', [180 40 200 35], 'string', 'Update downward ', 'fontsize', 13 ,...
    'BackgroundColor', [0.9 0 0], 'callback', @(Src,Ev)local_updatedown(hEdit));
set(figh,'toolbar', 'figure');

%=====================================================
function local_replot(hE, figh);
Thr = str2num(get(hE, 'string'));
Uidx = getGUIdata(gcbf, 'Uidx');
clf; drawnow;
JLspikeThr(Uidx, Thr, figh);

function local_rmcache(hE, figh);
Thr = str2num(get(hE, 'string'));
Uidx = getGUIdata(gcbf, 'Uidx');
JLspikeThr(Uidx, nan, 'rmcache');
close(gcbf);
drawnow;
JLspikeThr(Uidx, Thr);

function local_update(hE);
Thr = str2num(get(hE, 'string'));
Uidx = getGUIdata(gcbf, 'Uidx');
JL_APparam(Uidx, 'APthrSlope', Thr);
JL_APparam(Uidx, 'CutoutThrSlope', Thr); % also remove cache, so JLwaveforms will force new read action
JLwaveforms(Uidx,'rmcache');
close(gcbf);

function local_updatedown(hE);
Thr = str2num(get(hE, 'string'));
Uidx = getGUIdata(gcbf, 'Uidx');
DB = JLdbase;
db = JLdbase(Uidx);
DB = DB([DB.icell_run]==db.icell_run); % all recs of this cell ...
DB = DB([DB.MinutesSinceRecStart]>=db.MinutesSinceRecStart); % ... starting w current rec
UIDX = [DB.UniqueRecordingIndex];
for ii=1:numel(UIDX),
    JL_APparam(UIDX(ii), 'APthrSlope', Thr);
    JL_APparam(UIDX(ii), 'CutoutThrSlope', Thr);
    JLwaveforms(UIDX(ii),'rmcache'); % also remove cache, so JLwaveforms will force new read action
end
close(gcbf);









