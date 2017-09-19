function Y = JLbeatVar(S, doPlot);
% JLbeatVar - analyze variance of response to stimulus periods
%    JLbeatVar(S, doPlot);

if nargin<2,
    doPlot = nargout<1;
end

preDur = 500; % ms pre-stim by convention
postDur = 500; % ms post stim after which to analyze spont
tstartSpont = preDur+S.burstDur+postDur; % evaluate spont activity starting 500 ms after stim offset
Ttrans = S.Ttrans; % transition time after onset, excluded from analysis

Uidx = S.UniqueRecordingIndex;
%S = JLbeatList(S.ExpID, S.icell);
S = S([S.UniqueRecordingIndex]==Uidx);

[D, DS, L]=readTKABF(S);
dt = D.dt_ms;
S.dt = dt;
rec = D.AD(1).samples; % recording
NsamTot = D.NsamPerSweep/D.Nchan;
cutWin = [S.APwindow_start S.APwindow_end];
[rec, SPTraw, Tsnip, Snip, APslope] = APtruncate2(rec, S.CutoutThrSlope, dt, cutWin, 1); % last 1: replace APs by nans
%hist(rec, 100); return
%dplot(D.dt_ms, rec); return
clear D;

% spontaneous waveform: after beat
isamSpont = 1+round(tstartSpont/dt):NsamTot;
if numel(isamSpont)*dt<preDur,  % post-stim part is too short; take spont trace from start of recording 
    tstartSpont = 0;
    isamSpont = 1:round(preDur/dt);
end
Vspont = rec(isamSpont);
SpontMean = nanmean(Vspont);
SpontVar = nanvar(Vspont);

% analysis window, beat counting, etc
tstart = preDur + Ttrans; % ms start of analysis
tend = preDur + S.burstDur; % ms end of analysis
BeatCycle = 1e3/S.Fbeat; % ms beat duration
NbeatCycle = floor((tend-tstart)/BeatCycle); % # beats analyzed
NsamBeat = round(BeatCycle/dt); % # samples in one beat 
isam = round(tstart/dt) + (1:(NbeatCycle*NsamBeat)); % sample range
tend = dt*(isam(end)-1); % end of analysis revisited: integer # beat cycles
rec = rec(isam); % apply analysis window

% APs
SPT = SPTraw(betwixt(SPTraw, tstart,tend)); % select spikes occurring during analysis window
SPT = rem(SPT-preDur, BeatCycle);  %spike times re beat cycle start

% periodicity analysis: ipsi, contra, bin & spont
%function [VA, rec] = local_VarAna(rec, S, Cycle, preDur, sep);
Period1 = 1e3/S.Freq1; Period2 = 1e3/S.Freq2;
[VA_ipsi, rec_i] = local_VarAna(rec, S, Period1, Ttrans, 'ipsi'); 
[VA_contra, rec_c] = local_VarAna(rec, S, Period2, Ttrans, 'contra');
VA_bin = local_VarAna(rec, S, [BeatCycle Period1], Ttrans ,'bin');
VA_sil = local_VarAna(Vspont, S, Period1, preDur, 'spont');

VA_ipsi = AddPrefix(VA_ipsi, 'ipsi_');
VA_contra = AddPrefix(VA_contra, 'contra_');
VA_bin.PctgTimeBelowSpontVar = 100*nanmean(VA_bin.VarWave<SpontVar);
VA_bin = AddPrefix(VA_bin, 'bin_');

VA_sil = AddPrefix(VA_sil, 'sil_');

% ---assemble return arg Y---
Y = S;
Y = structJoin(Y, CollectInStruct(Period1, Period2));
sep = struct('x_____________________x123', '__________________');
Y = structJoin(Y, sep, CollectInStruct(preDur, Ttrans, tstart, tend, BeatCycle, NbeatCycle, sep, SPTraw, SPT));
Y.SpontVar = SpontVar;
Y.SpontMean = SpontMean;
sep = struct('x_____________________VVVVV', '__________________');
Y = structJoin(Y, sep, VA_ipsi, VA_contra, VA_bin, VA_sil);

[dt_pred, MeanPred, VarPred] = local_sumpred(Y);
Y.CorrMeanPred = corr(MeanPred, Y.bin_MeanWave);
Y.CorrVarPred = corr(VarPred, Y.bin_VarWave);

[Y.MaxCC_ipsi, ilag] = maxcorr(Y.ipsi_VarWave, Y.ipsi_MeanWave);
Y.tau_ipsi = ilag*Y.ipsi_dt;
[Y.MaxCC_contra, ilag] = maxcorr(Y.contra_MeanWave, Y.contra_VarWave);
Y.tau_contra = ilag*Y.contra_dt;

%==================================
%==================================

if doPlot,
    ThickGray = struct('color', 0.7*[1 1 1], 'linewidth', 3);
    %
    % ipsi and contra cycles superimposed & cyle-averaged variances
%     fh = figure;
%     set(fh,'units', 'normalized', 'position', [0.666 0.104 0.352 0.554]);
%     figname(fh, [mfilename 'ICvar']);
%     JLfigmenu(fh,S);
%     JLfigmenu(fh, S);
%     subplot(2,1,1); %-----ipsi traces--------
%     dplot(Y.ipsi_dt, rec_i);
%     xdplot(Y.ipsi_dt, Y.ipsi_MeanWave, 'k', 'linewidth',3);
%     title(S.TTT);
%     subplot(2,1,2);  %-----contra traces--------
%     dplot(Y.contra_dt, rec_c);
%     xdplot(Y.contra_dt, Y.contra_MeanWave, 'k', 'linewidth',3);
%     xlabel('Time (ms)');
    %
    % BEAT cycle mean & variance. Spont values as reference.
    fh2 = figure;
    set(fh2,'units', 'normalized', 'position', [0.0203 0.646 0.972 0.295]);
    figname(fh2, [mfilename 'Beatvar'])
    JLfigmenu(fh2,S);
    subplot(2,1,1); % variance
    dplot(Y.bin_dt, Y.bin_VarWave);
    xplot(xlim, SpontVar*[1 1], ThickGray);
    xdplot(Y.bin_dt, Y.bin_VarWave, 'k', 'linewidth', 2);
    xdplot(dt_pred, VarPred, 'm');
    title([S.TTT '  ' sprintf('binvar<spontvar %i%% of time', round(Y.bin_PctgTimeBelowSpontVar))]);
    subplot(2,1,2); % mean
    local_featherplot(Y.bin_dt, Y.bin_MeanWave, Y.bin_VarWave, [0 0.6 0], [0 0 0]);
    xdplot(dt_pred, MeanPred, 'm');
    fenceplot(SPT, min(ylim)+[0 0.25*diff(ylim)], 'color', [1 0.5 0.5]);
    xlabel('Time (ms)');
    TracePlotInterface(fh2);
    %
    % cycled-averaged mean vs ditto variance
    set(figure,'units', 'normalized', 'position', [0.291 0.219 0.348 0.357]);
    figname(gcf, [mfilename 'mean_vs_var']);
    plot(Y.bin_MeanWave, Y.bin_VarWave, '.', 'color', 0.9*[0.8 1 0.8], 'markersize', 4);
    xplot(Y.ipsi_MeanWave, Y.ipsi_VarWave, 'b.');
    xplot(Y.contra_MeanWave, Y.contra_VarWave, 'r.');
    xplot(Y.sil_MeanWave, Y.sil_VarWave, '.', 'color', 0.7*[1 1 1]);
    xplot(mean(Y.ipsi_MeanWave), mean(Y.ipsi_VarWave), 'b+', 'markersize', 10);
    xplot(mean(Y.contra_MeanWave), mean(Y.contra_VarWave), 'r+', 'markersize', 10);
    xplot(mean(Y.sil_MeanWave), mean(Y.sil_VarWave), '+', 'color', 0.7*[1 1 1], 'markersize', 10);
    xlabel('Mean (mV)');
    ylabel('Variance (mV^2)');
    ylim([0 max(ylim)]);
    %
    % cycled-averaged mean & variance vs time
    set(figure,'units', 'normalized', 'position', [0.0305 0.164 0.248 0.422]);
    figname(gcf, [mfilename 'meanvar_vs_time']);
    subplot(2,1,2); % mean
    local_featherplot(Y.ipsi_dt, Y.ipsi_MeanWave, Y.ipsi_VarWave, [0 0 1]);
    local_featherplot(Y.contra_dt, Y.contra_MeanWave, Y.contra_VarWave, [1 0 0]);
    xplot(xlim, [1 1]*Y.SpontMean, ThickGray);
    ylabel('cycle-av. waveform');
    xlabel('Time (ms)');
    %legend({'ipsi' 'contra'});
    subplot(2,1,1); % var
    dplot(Y.ipsi_dt, Y.ipsi_VarWave, 'b');
    xdplot(Y.contra_dt, Y.contra_VarWave, 'r');
    xplot(xlim, [1 1]*Y.SpontVar, ThickGray)
    xplot(xlim, mean(Y.ipsi_VarWave)*[1 1], 'b:')
    xplot(xlim, mean(Y.contra_VarWave)*[1 1], 'r:')
    ylabel('cycle-av. variance');
    ylim([0 max(ylim)]);
    title(sprintf('CC=[%1.2f, %1.2f] lag_{var}=[%i, %i] \\mus', ...
        Y.MaxCC_ipsi, Y.MaxCC_contra, round(Y.tau_ipsi*1e3), round(Y.tau_contra*1e3)));
    equalize_tagged_figs XY;
    %
    % sum prediction versus beat
    fh5 = figure;
    set(fh5,'units', 'normalized', 'position', [0.661 0.268 0.338 0.201]);
    figname(fh5, [mfilename 'beatrec_vs_sumpred']);
    subplot(1,2,1);
    plot(MeanPred, Y.bin_MeanWave, 'k');
    xlim([min([xlim ylim])   max([xlim ylim])]); ylim(xlim); axis square;
    xplot(xlim, ylim, 'k:');
    xplot(MeanPred, Y.bin_MeanWave, 'b');
    title(sprintf('rho = %1.2f', Y.CorrMeanPred));
    xlabel('Measured V_m (mV)');
    ylabel('Predicted V_m (mV)');
    subplot(1,2,2);
    plot(VarPred, Y.bin_VarWave, 'k');
    xlim([min([xlim ylim])   max([xlim ylim])]); ylim(xlim); axis square;
    xplot(xlim, ylim, 'k:');
    xplot(VarPred, Y.bin_VarWave, 'b');
    title(sprintf('rho = %1.2f', Y.CorrVarPred));
    xlabel('Measured V_m (mV)');
    ylabel('Predicted V_m (mV)');
    %
    figure(fh2);
end


%===================================
function [VA, rec] = local_VarAna(rec, S, Cycle, Ttrans, sep);
%RejectAPcycles = arginDefaults('RejectAPcycles',0);
if numel(Cycle)==1,
    smallCycle=0;
else
    [Cycle, smallCycle] = DealElements(Cycle);
end
Nsam = numel(rec);
dt = S.dt;
Dur = Nsam*dt;
Ncycle = floor(Dur/Cycle); % # cycles considered
Tend = Ncycle*Cycle; % time of last point of rec used for loop mean
Nsam = round(Tend/dt);
rec = rec(1:Nsam); % truncate tail
NsamCycle = round(Nsam*Cycle/Dur); % target #samples per cycle
NewNsam = Ncycle*NsamCycle; % new # samples, for which cycle has integer # points
warning('off', 'MATLAB:interp1:NaNinY');
rec = interp1((1:Nsam).',rec,linspace(1,Tend/dt,NewNsam).'); % primitive resample by linear interp
warning('on', 'MATLAB:interp1:NaNinY');
dt = dt*Nsam/NewNsam; % new sample rate
if isequal(0, smallCycle), % no double periodicity
    rec = reshape(rec, [NsamCycle Ncycle]);
else, % beat cycle may be incommensurate w tone cycle. Correct for this mismatch.
    rec = local_reshape(rec, NsamCycle, Ncycle, dt, smallCycle);
end
NsamShift = round(Ttrans/dt); % shift rec so that t=0 corresponds to stim onset
if isequal(0, smallCycle), % no double periodicity
    rec = circshift(rec,NsamShift); % a delay that matches the delay btw stim onset and onset of ana window
else,  % beat cycle may be incommensurate w tone cycle. Circshift using a shortened beat cycle
    rec = local_circshift(rec, NsamShift, dt, smallCycle);
end
AC = shufcorr(rec(:,all(~isnan(rec)))); % only consider columns free of NaNs
MeanWave = nanmean(rec,2);
VarWave = nanvar(rec,0,2);
totVar = nanvar(rec(:));
resVar = mean(VarWave);
varAcc = 100*((totVar-resVar)/totVar);
x______________ = '____________';
VA = CollectInStruct(x______________, dt, AC, MeanWave, VarWave, totVar, resVar, varAcc);

function rec = local_rejectCycles(rec, S, SPTraw, tstart, dt);
% reject cycles (columns) containing an AP
qq = 0*rec(:);
TT = timeaxis(qq, dt, tstart);
APstart = SPTraw + S.APwindow_start;
APend = SPTraw + S.APwindow_end;
N_AP = numel(APstart);
for iAP=1:N_AP,
    %qq(betwixt(TT,APstart(iAP),APend(iAP))) = 1;
    qq(TT>APstart(iAP) & TT<APend(iAP)) = 1;
end
%f5; dplot(dt, rec(:)); xdplot(dt, qq, 'r'); tracePlotInterface(gcf);
qq = reshape(qq,size(rec)); % each column corresponds to a cycle of rec
iok = ~any(qq); % select columns not containing an AP
%f6; dplot(dt, rec(:,iok)); xdplot(dt, rec(:,~iok), 'r'); 
rec = rec(:,iok);


function local_featherplot(dt, Mn, Vr, CLR, lineCLR);
lineCLR = arginDefaults('lineCLR', CLR);
Tm = Xaxis(Mn, dt);
St = sqrt(Vr);
Pstd = [Mn(:)+St(:) ; flipud(Mn(:)-St(:))];
Ptime = [Tm(:) ; flipud(Tm(:))];
gCLR = 0.85*(CLR+1)/max(CLR+1);
patch(Ptime, Pstd, gCLR, 'FaceAlpha', 0.5, 'edgecolor', gCLR);
xdplot(dt, Mn, 'linewidth', 2, 'color', lineCLR);


function  [dt, MeanPred, VarPred] = local_sumpred(Y);
% predict binbeat wave from monaural contribs
MeanPred = local_repcycle(Y.bin_dt, Y.bin_MeanWave, Y.ipsi_dt, Y.ipsi_MeanWave) ...
    + local_repcycle(Y.bin_dt, Y.bin_MeanWave, Y.contra_dt, Y.contra_MeanWave);
MeanPred = MeanPred - 0.5*mean(MeanPred);
VarPred = local_repcycle(Y.bin_dt, Y.bin_VarWave, Y.ipsi_dt, Y.ipsi_VarWave) ...
    + local_repcycle(Y.bin_dt, Y.bin_VarWave, Y.contra_dt, Y.contra_VarWave);
VarPred = VarPred - 0.5*mean(VarPred);
dt = Y.bin_dt;

function Ymon = local_repcycle(dt_bin, Ybin, dt_mon, Ymon);
BinDur = dt_bin*numel(Ybin);
MonDur = dt_mon*numel(Ymon);
NmonCycle = ceil(BinDur/MonDur)+1; % one cycle margin to avoid interp1 nans
Ymon = Ymon*ones(1,NmonCycle); % faster than repmat ;)
monTime = Xaxis(Ymon(:), dt_mon);
binTime = Xaxis(Ybin(:), dt_bin);
Ymon = interp1(monTime, Ymon(:), binTime);

function  rr = local_reshape(rec, NsamCycle, Ncycle, dt, SmallCycle);
%NsamCycle, Ncycle, dt, SmallCycle
tstart = dt*NsamCycle*(0:Ncycle-1); % ideal starting points for loop mean. ..
% ... The only problem is that some tstart values may not be integer
% multiples of SmallCycle. So round them to the nearest smaller multiples.
tstart = SmallCycle*floor(tstart/SmallCycle); % round down to avoid out-of-range errors
istart = 1+round(tstart/dt); % time -> sample index
iend = istart + NsamCycle -1;
rr = zeros(NsamCycle,Ncycle);
for icycle=1:Ncycle,
    rr(:,icycle) = rec(istart(icycle):iend(icycle));
end

function  rec = local_circshift(rec, NsamShift, dt, smallCycle);
% like circshift but accounting for non-integer # stim cycles in beat cycle
Nsam = size(rec,1);
Dur = dt*Nsam;
Ncycle = floor(Dur/smallCycle);
NsamShort = round(Ncycle*smallCycle/dt);
rec = rec(1:NsamShort,:);
rec = circshift(rec, NsamShift);
% add the lost samples by repeating the first fraction
rec = [rec; rec(1:(Nsam-NsamShort),:)];
















