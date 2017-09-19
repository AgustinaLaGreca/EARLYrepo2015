function [S, rawrec] = JLgetRec(Uidx, APthr, APwin);
% JLgetRec - get JLbeat recording w spikes & 50 Hz noise removed.
%   [S, Rawrec] = JLgetRec(Uidx) uses default AP param values as given by
%   JLdatastruct(uidx). S is a struct containing stimulus parameters,
%   processed waveforms (APs removed), spike waveforms and spike times, etc. 
%   Rawrec is the raw waveform with 50-Hz filtered out - no further processing.
%
%   [S, rawrec] = JLgetRec(Uidx, APthr, APwin) 
%   uses the specified AP param values APthr and APwin. 
%     APthr is the negative thresholds APthrSlope in V/s 
%       (typically about -10 V/s) or a 2-component array
%       [APthrSlope, CutoutThrSlope]. Specify APthr==-inf to leave APs alone.
%     APwin is 2-comp array with onset and offset of APs re their peak.
%       A typical value is [-1 1.4] ms. Details on AP detection are given in
%       APtruncate2.
% 
%   Caching is used only if rawrec is not requested and no AP params are
%   specified.
%
%   See also JLviewrec, APtruncate2.

[APthr, APwin] = arginDefaults('APthr/APwin', []); % default: take from JLdatastruct info (struct S defined below)

CFN = fullfile(processed_datadir, 'JL', 'qraw', ['Rec__' num2str(Uidx) '.cache']);
DoCache = nargin<2; % default AP params -> use caching

S = JLdatastruct(Uidx);
APP = JL_APparam(Uidx);
preDur = 500; % ms pre-stim by convention
postDur = 450; % ms post stim after which to analyze spont
Ttrans = 500; % ms post-onset transition
Segmentation = CollectInStruct(preDur, postDur, Ttrans);

if isempty(APthr), % use values in APP
    [APthrSlope, CutoutThrSlope] = deal(APP.APthrSlope, APP.CutoutThrSlope);
else, % overrule the values in APP
    [APthrSlope, CutoutThrSlope] = deal(APthr(1), APthr(end));
end
if isempty(APwin), 
    [APwindow_start APwindow_end] = deal(APP.APwindow_start, APP.APwindow_end);
else, % overrule the values in APP
    [APwindow_start APwindow_end] = deal(APwin(1), APwin(2));
end
% make sure AP param values in S match those used here
[S.APthrSlope, S.CutoutThrSlope, S.APwindow_start, S.APwindow_end] ...
    = deal(APthrSlope, CutoutThrSlope, APwindow_start, APwindow_end);
% ditto Segment
[S.preDur, S.postDur, S.Ttrans] ...
    = deal(Segmentation.preDur, Segmentation.postDur, Segmentation.Ttrans);

% Only now that we fixed AP params & Segmentation, try to retrieve any cached value
if DoCache,
    [s, CFN, CP] = getcache(CFN, {APP Segmentation});
    if ~isempty(s) && nargout<2,
        S = s;
        return;
    end
else, % we're probably testing alternativeAP thresholds. Don't cache these attempts.
end


SpontLoop = 10*sqrt(2); % ms interval over which to evaluate spont act
%APthrSlope, CutoutThrSlope, APwindow_start, APwindow_end

[D, DS, L]=readTKABF(S);
dt = D.dt_ms;
rawrec = D.AD(1).samples; % entire recording, including pre- and post-stim parts
rawrec = local_reject50Hz(rawrec, dt);
%rec = D.AD(3).samples + D.AD(4).samples;  %DEBUG  stimuli of the two ears summed
NsamTot = numel(rawrec);
cutWin = [APP.APwindow_start APP.APwindow_end];
[rec, SPTraw, Tsnip, Snip, APslope] = ...
    APtruncate2(rawrec, [APthrSlope, CutoutThrSlope], dt, cutWin, 1); % last 1: replace APs by nans
clear D;

% spontaneous waveform: after beat
tstartSpont1 = 0; % ms start of first spont interval
isamSpont1 = 1:round(preDur/dt);
tstartSpont2 = preDur+S.burstDur+postDur; % evaluate spont activity starting 500 ms after stim offset
isamSpont2 = 1+round(tstartSpont2/dt):NsamTot;
NsamLoop = round(SpontLoop/dt);
[SpontDur1, SpontMean1, SpontVar1, FastSpontVar1] = local_spont(rec(isamSpont1), NsamLoop, dt);
[SpontDur2, SpontMean2, SpontVar2, FastSpontVar2] = local_spont(rec(isamSpont2), NsamLoop, dt);


% analysis window, beat counting, etc
tstart = preDur + Ttrans; % ms start of analysis
tend = preDur + S.burstDur; % ms end of analysis
BeatCycle = 1e3/S.Fbeat; % ms beat duration
NbeatCycle = floor((tend-tstart)/BeatCycle); % # beats analyzed
NsamBeat = round(BeatCycle/dt); % # samples in one beat
isam = round(tstart/dt) + (1:(NbeatCycle*NsamBeat)); % sample range
rec = rec(isam); % apply analysis window

recDur = numel(rec)*dt; % duration of rec

% ====spike times====
SPT = SPTraw-preDur; % spike times expressed re stim onset
SPT = SPT(betwixt(SPT, Ttrans, Ttrans+recDur)); % only consider spikes occurring within analysis window

S.Spikes___________________ = '_____________________';
S.SPT = SPT;
S.SPTraw = SPTraw;
S.Tsnip = Tsnip;
S.Snip = Snip;

qq.Spont____________________ = '_____________________';
qq = structJoin(qq, CollectInStruct(tstartSpont1, tstartSpont2, SpontLoop, ...
    SpontDur1, SpontMean1, SpontVar1, FastSpontVar1, ...
    SpontDur2, SpontMean2, SpontVar2, FastSpontVar2));
S = structJoin(S, qq);

qqq.Recording________________ = '_____________________';
qqq = structJoin(qqq, CollectInStruct(dt, preDur, postDur, recDur, rec));
S = structJoin(S, qqq);

if DoCache, 
    putcache(CFN, 5e3, CP, S); 
end


%============================================================
function [Dur, Mn, Vr, FastVr] = local_spont(rec, NsamLoop, dt);
Dur = numel(rec)*dt;
Mn = nanmean(rec(:));
Vr = nanvar(rec(:));
N = numel(rec)-rem(numel(rec), NsamLoop);
rec = reshape(rec(1:N), NsamLoop, N/NsamLoop);
FastVr = nanmean(nanvar(rec,1));

function y = local_reject50Hz(y, dt);
freq = 50; % Hz
Z = cosfit(dt, y-mean(y), freq);
RunPhase = 2*pi*freq*Xaxis(y,dt/1000) + angle(Z);
y = y - abs(Z)*cos(RunPhase);


