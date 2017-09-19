function [T, W, RecMat] = JLviewrec(Uidx, doPlot);
% JLviewrec - view raw data with basic feature detection
%   JLviewrec(Uidx) plots
%   [T, W] = JLviewrec(Uidx, doPlot) returns stats T and waveforms W.
%   Waveforms are read from raw data and subjected to spike thresholding.
%
%   See also JLgetrec.

doPlot = arginDefaults('doPlot',nargout<1);

[S, rawrec] = JLgetRec(Uidx);
dt = S.dt;
isam = @(t)1+round(t/dt);
titv = @(T)dt*(isam(T(1)):isam(T(2))).'-dt;
BeatDur = 1e3/S.Fbeat; % ms beat cycle duration
NbeatCycle = floor((S.burstDur-S.Ttrans)/BeatDur); % # beat cycles used in loopmean
TotDur = dt*numel(rawrec);

% slow trend
DTslow = BeatDur; % ms sampe period for trend 
Ndec = round(DTslow/dt); % decimation factor
GMean = mean(rawrec(:)); % grand mean
recSlow = resample(rawrec-GMean, 1, Ndec);
recSlow = GMean+resample(recSlow, Ndec, 1);
Nhalf = round(Ndec/10); % decimation factor for somewhat slow trends
recHalfSlow = GMean+resample(rawrec-GMean, 1, Nhalf);
dt_half = dt*Nhalf;

% bookkeeping & Landmarks 
NsamLoop = numel(S.rec);
NsamTot = numel(rawrec);
tonset = S.preDur; % ms stim onset
tstart = S.preDur + S.Ttrans; % ms start of loop mean analysis 
tend = tstart + NbeatCycle*BeatDur; % end of loopmean analysis
toffset = S.preDur + S.burstDur; % ms stim offset
tspont2 = toffset + S.postDur; % ms onset of post-stim spont interval
AnaDur = tend-tstart; % ms duration of analysis 

% global means; needed now for recovering DC lost by detrending below
rr = cutFromWaveform(dt, rawrec, S.SPTraw, [S.APwindow_start S.APwindow_end], 1); % last 1: turn spikes into nans
MeanSpont1 = nanmean(rr(1:isam(tonset))); % var of 1st spont itv
MeanOnset = nanmean(rr(isam(tonset):isam(tstart)-1));
MeanDriv = nanmean(rr(isam(tstart):isam(tend)-1));
MeanTail = nanmean(rr(isam(tend):isam(toffset)-1));
MeanOffset = nanmean(rr(isam(toffset):isam(tspont2)-1));
MeanSpont2 = nanmean(rr(isam(tspont2):end));
Spont2Dur = TotDur-tspont2;

% Waveforms. The detrended waveform during the stimulus (after adapting
% for S.Ttrans ms) serves as the basis for the various cycle means.
recSlow_anawin = recSlow(isam(titv([tstart tend-dt])));
stimrec_detrended = MeanDriv + S.rec-recSlow_anawin; % the detrended waveform during stimulus. DC re-inserted
[dt_BeatMean, BeatMeanrec, BeatVar, RecMatBeat] = fracLoopMean(dt, stimrec_detrended, BeatDur, tstart-tonset); % beat-cycle loop mean of analysis window
isamloop = @(t)1+round(t/dt_BeatMean); 
Tstim = titv([tonset toffset]); % time axis of stimulus 
TreBeat = mod(Tstim-tstart, BeatDur); % time re beat cycle 
IpsiCycle = 1e3/S.Freq1; % ms ipsi stim cyle
[dt_IpsiMean, IpsiMeanrec, IpsiVar, RecMatIpsi] = fracLoopMean(dt, stimrec_detrended, IpsiCycle, tstart-tonset); % ipsi-cycle loop mean of analysis window
ContraCycle = 1e3/S.Freq2; % ms contra stim cyle
[dt_ContraMean, ContraMeanrec, ContraVar, RecMatContra] = fracLoopMean(dt, stimrec_detrended, ContraCycle, tstart-tonset); % ipsi-cycle loop mean of analysis window
% peak values of cycle means
MaxMeanIpsiCycle = max(IpsiMeanrec);
MinMeanIpsiCycle = min(IpsiMeanrec);
MaxVarIpsiCycle = max(IpsiVar);
MinVarIpsiCycle = min(IpsiVar);
MaxMeanContraCycle = max(ContraMeanrec);
MinMeanContraCycle = min(ContraMeanrec);
MaxVarContraCycle = max(ContraVar);
MinVarContraCycle = min(ContraVar);
P2PMeanIpsiCycle = MaxMeanIpsiCycle - MinMeanIpsiCycle;
P2PVarIpsiCycle = MaxVarIpsiCycle - MinVarIpsiCycle;
P2PMeanContraCycle = MaxMeanContraCycle - MinMeanContraCycle;
P2PVarContraCycle = MaxVarContraCycle - MinVarContraCycle;
% spike waveforms
Nsp = size(S.SPTraw,1);
Tsp = nan; Wsp = nan;
for ii=1:Nsp,
    SpikeItv = S.SPTraw(ii)+[S.APwindow_start S.APwindow_end];
    tim = titv(SpikeItv); % time axis of spike
    tim = tim(betwixt(tim,0,TotDur-dt));
    Tsp = [Tsp; tim; nan]; % time axis
    Wsp = [Wsp; ; rawrec(isam(tim)); nan]; % waveform
end

% var: detrend first to eliminate contrib of slow trends
rr = cutFromWaveform(dt, rawrec-recSlow+GMean, S.SPTraw, [S.APwindow_start S.APwindow_end], 1); % last 1: turn spikes into nans
VarSpont1 = nanvar(rr(1:isam(tonset))); % var of 1st spont itv
VarOnset = nanvar(rr(isam(tonset):isam(tstart)-1));
VarDriv = nanvar(rr(isam(tstart):isam(tend)-1));
VarTail = nanvar(rr(isam(tend):isam(toffset)-1));
VarOffset = nanvar(rr(isam(toffset):isam(tspont2)-1));
VarSpont2 = nanvar(rr(isam(tspont2):end));
% slow variance (drift)
rr = cutFromWaveform(dt, recSlow, S.SPTraw, [S.APwindow_start S.APwindow_end], 1); % last 1: turn spikes into nans
VarDriftSpont1 = nanvar(rr(isam(tonset):isam(tstart)-1));
VarDriftOnset = nanvar(rr(isam(tonset):isam(tstart)-1));
VarDriftDriv = nanvar(rr(isam(tstart):isam(tend)-1));
VarDriftTail = nanvar(rr(isam(tend):isam(toffset)-1));
VarDriftOffset = nanvar(rr(isam(toffset):isam(tspont2)-1));
VarDriftSpont2 = nanvar(rr(isam(tspont2):end));
% spike waveforms revisited; use detrended waveforms, but insert mean
% true potentials afterwards, while distinguishing driven & spont portions
% of the recording.
[rr, Tsnip, Snip] = cutFromWaveform(dt, rawrec-recSlow, S.SPTraw, [S.APwindow_start S.APwindow_end], 1); % last 1: turn spikes into nans
SPTraw = S.SPTraw(:).';
APinSpont1 = betwixt(SPTraw, 0, tonset);
APinOnset = betwixt(SPTraw, tonset, tstart);
APinStim = betwixt(SPTraw, tstart, tend);
APinTail = betwixt(SPTraw, tend, toffset);
APinOffset = betwixt(SPTraw, toffset, tspont2);
APinSpont2 = betwixt(SPTraw, tspont2, inf);
Snip(:, APinSpont1) = Snip(:, APinSpont1) + MeanSpont1;
Snip(:, APinOnset) = Snip(:, APinOnset) + MeanOnset;
Snip(:, APinStim) = Snip(:, APinStim) + MeanDriv;
Snip(:, APinTail) = Snip(:, APinTail) + MeanTail;
Snip(:, APinOffset) = Snip(:, APinOffset) + MeanOffset;
Snip(:, APinSpont2) = Snip(:, APinSpont2) + MeanSpont2;
NspikeTotal = numel(APinStim);
NspikeInSpont1 = sum(APinSpont1);
NspikeInOnset = sum(APinOnset);
NspikeInStim = sum(APinStim);
NspikeInTail = sum(APinTail);
NspikeInOffset = sum(APinOffset);
NspikeInSpont2 = sum(APinSpont2);
Spont1Dur = S.preDur;
SpontRate1 = 1e3*NspikeInSpont1/Spont1Dur; % sp/s
SpontRate2 = 1e3*NspikeInSpont2/Spont2Dur; % sp/s
DrivenRate = 1e3*NspikeInStim/AnaDur; % sp/s

% return args
[t_anaStart, t_anaEnd] = deal(tstart, tend);
[t_stimOnset, t_stimOffset] = deal(tonset, toffset);
t_spont2Start = tspont2;
% residual vars
IpsiResVar = mean(IpsiVar);
ContraResVar = mean(ContraVar);
BeatResVar = mean(BeatVar);
UniqueRecordingIndex = S.UniqueRecordingIndex;
T = CollectInStruct(UniqueRecordingIndex, ...
    '-dur', TotDur, AnaDur, Spont1Dur, Spont2Dur, ...
    t_stimOnset, t_anaStart, t_anaEnd, t_stimOffset, t_spont2Start, ...
    '-mean',        MeanSpont1, MeanOnset, MeanDriv, MeanTail, MeanOffset, MeanSpont2, ...
    '-peaks_troughs', MaxMeanIpsiCycle, MinMeanIpsiCycle, MaxVarIpsiCycle, MinVarIpsiCycle, ... 
    MaxMeanContraCycle, MinMeanContraCycle, MaxVarContraCycle, MinVarContraCycle, ...
    P2PMeanIpsiCycle, P2PVarIpsiCycle, P2PMeanContraCycle, P2PVarContraCycle, ...
    '-resvar', IpsiResVar, ContraResVar, BeatResVar, ...
    '-detrended_var', VarSpont1, VarOnset, VarDriv, VarTail, VarOffset, VarSpont2, ...
    '-slow_var', VarDriftSpont1, VarDriftOnset ,VarDriftDriv, VarDriftTail, VarDriftOffset, VarDriftSpont2, ...
    '-spikestats', ...
    NspikeTotal, NspikeInSpont1, NspikeInOnset, NspikeInStim, NspikeInTail, NspikeInOffset, NspikeInSpont2, ...
    DrivenRate, SpontRate1, SpontRate2);
W = CollectInStruct(UniqueRecordingIndex, ...
    '-dur', TotDur, AnaDur, Spont1Dur, Spont2Dur, ...
    t_anaStart, t_anaEnd, t_stimOnset, t_stimOffset, t_spont2Start, ...
    '-waveforms', dt, BeatDur, stimrec_detrended, ...
    dt_BeatMean, BeatMeanrec, BeatVar, dt_IpsiMean, IpsiMeanrec, IpsiVar, ...
    dt_ContraMean, ContraMeanrec, ContraVar, ...
    '-spikes', SPTraw, Tsnip, Snip, APinSpont1, APinOnset, APinStim, APinTail, APinOffset, APinSpont2);
if nargout>2,
    [dt_Beat, dt_Ipsi, dt_Contra] = deal(dt_BeatMean, dt_IpsiMean, dt_ContraMean);
    RecMat = CollectInStruct(UniqueRecordingIndex, ...
        dt_Beat, RecMatBeat, dt_Ipsi, RecMatIpsi, dt_Contra, RecMatContra);
end

if doPlot, % ==========plot=============
    set(gcf,'units', 'normalized', 'position', [0.0117 0.499 0.95 0.41]);
    hraw = xdplot(S.dt, rawrec, 'color', 0.7*[1 1 1], 'linewidth',2); % raw waveform.
    hdetr = xdplot(S.dt, rawrec-recSlow+GMean, 'color', 'k', 'linewidth',1.5); % detrended waveform.
    hloop = xplot(Tstim, BeatMeanrec(isamloop(TreBeat)), 'color', [0.1 0.7 0.1], 'linewidth', 2);
    htrend = xdplot(dt, recSlow, 'color', [0.3 0.3 1]);
    htrend2 = xdplot(dt_half, recHalfSlow, 'color', 0.5*[0.3 0.3 1]);
    hspike = xplot(Tsp, Wsp, 'r'); % spike waveforms
    xplot(S.SPTraw, rawrec(isam(S.SPTraw)), 'x', 'color', [0.8 0 0], 'markersize', 8); % spike times
    % landmarks
    ylim(ylim); % fix y limits
    fenceplot([tstart tend tspont2], ylim, 'color', 'm', 'linewidth', 2);
    fenceplot([tonset toffset], ylim, 'color', 'g', 'linewidth', 1.5);
    % rec info
    title(['exp ' num2str(S.iexp) '  cell ' num2str(S.icell) '   ' S.chan ' ' num2str(S.freq) ' Hz  '  num2str(S.SPL) ' dB'  ]);
    hl = legend([hraw hdetr htrend htrend2 hloop hspike], {'raw' 'detrended' 'trend' 'fast trend' 'Bt cyc mean' 'spikes'}, ...
        'location', 'SouthWest');
    LegPos = get(hl, 'Position'); LegPos([1 2]) = [0.001 0.03];
    set(hl, 'Position', LegPos);
    xlabel('Time (ms)', 'fontsize', 14);
    ylabel('Potential (mV)', 'fontsize', 14);
    TracePlotInterface(gcf);
end





