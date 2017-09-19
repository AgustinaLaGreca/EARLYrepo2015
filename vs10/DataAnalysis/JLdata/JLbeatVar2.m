function Y = JLbeatVar2(Uidx, doPlot);
% JLbeatVar2 - analyze variance of response to stimulus periods
%    JLbeatVar2(Uidx, doPlot);

[doPlot] = arginDefaults('doPlot', nargout<1);

%S = JLdatastruct(Uidx); % standardized recording info
D = JLdbase(Uidx); % general info on the recording & stimulus
S = JLdatastruct(Uidx); % specific info on stimulus etc
W = JLwaveforms(Uidx); % waveforms and duration params
T = JLcycleStats(Uidx); % cycle-based means & variances

dt = W.dt;

% ====some extra analyses=====
% lin predictions and their corr with real rec
MeanPred = localPred(dt, W.BeatMeanrec, W.dt_IpsiMean, W.IpsiMeanrec, W.dt_ContraMean, W.ContraMeanrec);
VarPred = localPred(dt, W.BeatVar, W.dt_IpsiMean, W.IpsiVar, W.dt_ContraMean, W.ContraVar);
%dsize(MeanPred, W.BeatMeanrec);
CorrMeanPred = corr(MeanPred, W.BeatMeanrec);
CorrVarPred = corr(VarPred, W.BeatVar);

% spike times re begin of beat cycle
SPTstim = W.SPTraw(W.APinStim);
SPTstim = mod(SPTstim - W.t_anaStart, W.BeatDur); 
% correlation & lag btw cycle variance and cycle mean
[MaxCC_ipsi, ilag] = maxcorr(W.IpsiVar, W.IpsiMeanrec);
tau_ipsi = ilag*W.dt_IpsiMean;
[MaxCC_contra, ilag] = maxcorr(W.ContraVar, W.ContraMeanrec);
tau_contra = ilag*W.dt_ContraMean;

% CorrMeanPred,  MaxCC_contra,  Uidx,        tau_contra    
% CorrVarPred,  MaxCC_ipsi, SPTstim,  VarPred, tau_ipsi      
% MeanPred, 
UniqueRecordingIndex = Uidx;
Y = CollectInStruct(UniqueRecordingIndex, '-cycleCorr', ...
    CorrMeanPred, CorrVarPred, MaxCC_ipsi, tau_ipsi, MaxCC_contra, tau_contra);
% primitive caching derived from JL_NNTP
local_cache(Y);
%==================================
%==================================

if doPlot,
    TT = JLtitle(Uidx);
    ThickGray = struct('color', 0.7*[1 1 1], 'linewidth', 3);
    %
    % ==========BEAT cycle mean & variance. Spont values as reference.
    fh2 = figure;
    set(fh2,'units', 'normalized', 'position', [0.0203 0.646 0.972 0.295]);
    figname(fh2, [mfilename 'Beatvar'])
    JLfigmenu(fh2,S);
    subplot(2,1,1); % variance
    dplot(W.dt_BeatMean, W.BeatVar);
    local_arrow(T.VarSpont1, '', 'left');
    local_arrow(T.VarSpont1, 'SPONT1', 'right');
    local_arrow(T.IpsiResVar, '', 'left', 'color', 'b');
    local_arrow(T.IpsiResVar, 'IPSI', 'right', 'color', 'b');
    local_arrow(T.ContraResVar, '', 'left', 'color', 'r');
    local_arrow(T.ContraResVar, 'CONTRA', 'right', 'color', 'r');
    %xplot(xlim, T.VarSpont1*[1 1], ThickGray);
    xdplot(W.dt_BeatMean, W.BeatVar, 'k', 'linewidth', 2);
    xdplot(W.dt_BeatMean, VarPred, 'm');
    ylabel('Var (mV^2)');
    title(TT);
    %title([S.TTT '  ' sprintf('binvar<spontvar %i%% of time', round(Y.bin_PctgTimeBelowSpontVar))]);
    subplot(2,1,2); % mean
    MM = local_featherplot(W.dt_BeatMean, W.BeatMeanrec, W.BeatVar, [0 0.6 0], [0 0 0]);
    ylim(mean(MM)+0.75*[-1 1]*diff(MM));
    xdplot(W.dt_BeatMean, MeanPred, 'm');
    fenceplot(SPTstim, min(ylim)+[0 0.25*diff(ylim)], 'color', [1 0.5 0.5]);
    xlabel('Time (ms)');
    ylabel('Avg potential (mV)');
    TracePlotInterface(fh2);
    %
    % ===========cycled-averaged mean vs ditto variance
    set(figure,'units', 'normalized', 'position', [0.291 0.219 0.348 0.357]);
    figname(gcf, [mfilename 'mean_vs_var']);
    plot(W.BeatMeanrec, W.BeatVar, '.', 'color', 0.9*[0.8 1 0.8], 'markersize', 4);
    xplot(W.IpsiMeanrec, W.IpsiVar, 'b.');
    xplot(W.ContraMeanrec, W.ContraVar, 'r.');
    %xplot(Y.sil_MeanWave, Y.sil_VarWave, '.', 'color', 0.7*[1 1 1]);
    xplot(mean(W.IpsiMeanrec), mean(W.IpsiVar), 'b+', 'markersize', 10);
    xplot(mean(W.ContraMeanrec), mean(W.ContraVar), 'r+', 'markersize', 10);
    xplot(T.MeanSpont1, T.VarSpont1, '+', 'color', 0.7*[1 1 1], 'markersize', 10);
    xlabel('Mean (mV)');
    ylabel('Variance (mV^2)');
    ylim([0 max(ylim)]);
    text(0.1, 0.9, sprintf('CC=[%1.2f, %1.2f] lag_{var}=[%i, %i] \\mus', ...
        MaxCC_ipsi, MaxCC_contra, round(tau_ipsi*1e3), round(tau_contra*1e3)), 'units', 'normalized');
    %
    % ===========cycled-averaged mean & variance vs time (monaural)
    set(figure,'units', 'normalized', 'position', [0.0305 0.164 0.248 0.422]);
    figname(gcf, [mfilename 'meanvar_vs_time']);
    subplot(2,1,2); % mean
    local_featherplot(W.dt_IpsiMean, W.IpsiMeanrec , W.IpsiVar, [0 0 1]);
    local_featherplot(W.dt_ContraMean, W.ContraMeanrec , W.ContraVar, [1 0 0]);
    xplot(xlim, [1 1]*T.MeanSpont1, ThickGray);
    ylabel('cycle-av. waveform');
    xlabel('Time (ms)');
    %legend({'ipsi' 'contra'});
    subplot(2,1,1); % var
    dplot(W.dt_IpsiMean, W.IpsiVar, 'b');
    xdplot(W.dt_ContraMean, W.ContraVar, 'r');
    xplot(xlim, [1 1]*T.VarSpont1, ThickGray)
    xplot(xlim, mean(W.IpsiVar)*[1 1], 'b:')
    xplot(xlim, mean(W.ContraVar)*[1 1], 'r:')
    ylabel('cycle-av. variance');
    title(TT);
    ylim([0 max(ylim)]);
    equalize_tagged_figs XY;
    %
    % ================sum prediction versus beat
    fh5 = figure;
    set(fh5,'units', 'normalized', 'position', [0.661 0.268 0.338 0.201]);
    figname(fh5, [mfilename 'beatrec_vs_sumpred']);
    subplot(1,2,1);
    plot(MeanPred, W.BeatMeanrec, 'k');
    xlim([min([xlim ylim])   max([xlim ylim])]); ylim(xlim); axis square;
    xplot(xlim, ylim, 'k:');
    xplot(MeanPred, W.BeatMeanrec, 'b');
    title(sprintf('rho = %1.2f', CorrMeanPred));
    xlabel('Measured V_m (mV)');
    ylabel('Predicted V_m (mV)');
    subplot(1,2,2);
    plot(VarPred, W.BeatVar, 'k');
    xlim([min([xlim ylim])   max([xlim ylim])]); ylim(xlim); axis square;
    xplot(xlim, ylim, 'k:');
    xplot(VarPred, W.BeatVar, 'b');
    title(sprintf('rho = %1.2f', CorrVarPred));
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


function MM = local_featherplot(dt, Mn, Vr, CLR, lineCLR);
lineCLR = arginDefaults('lineCLR', CLR);
Tm = Xaxis(Mn, dt);
St = sqrt(Vr);
Pstd = [Mn(:)+St(:) ; flipud(Mn(:)-St(:))];
Ptime = [Tm(:) ; flipud(Tm(:))];
gCLR = 0.85*(CLR+1)/max(CLR+1);
patch(Ptime, Pstd, gCLR, 'FaceAlpha', 0.5, 'edgecolor', gCLR);
xdplot(dt, Mn, 'linewidth', 2, 'color', lineCLR);
[MM(1), MM(2)] = minmax(Pstd);


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


function Pred = localPred(dt, wB, dtI, wI, dtC, wC);
T = Xaxis(wB,dt); % time
BeatDur = numel(wB)*dt;
wI = wI - mean(wI)/2;
wC = wC - mean(wC)/2;
% repeat ipsi cycles & resample to match beat waveform
cycI = dtI*numel(wI); % ipsi cyle
NcycleI = ceil(BeatDur/cycI);
wI = repmat(wI(:), [NcycleI, 1]); % repeated wI having at least beat cycle  dur
timeI = Xaxis(wI, dtI);
wI = interp1(timeI, wI, T); % resample wI using lin interp
% repeat contra cycles & resample to match beat waveform
cycC = dtC*numel(wC); % contra cyle
NcycleC = ceil(BeatDur/cycC);
wC = repmat(wC(:), [NcycleC, 1]); % repeated wC having at least beat cycle  dur
timeC = Xaxis(wC, dtC);
wC = interp1(timeC, wC, T); % resample wC using lin interp
Pred = reshape(wI + wC, size(wB));

function local_arrow(Y, Str, loc, varargin);
Ynorm = interp1(ylim, [0 1], Y);
if isequal('left', loc),
    har = text(0, Ynorm, [Str '\rightarrow'], 'units', 'normalized', 'fontsize', 10, ...
        'fontweight', 'bold', 'horizontalalign', 'right', 'verticalalign', 'middle', varargin{:});
    %     XX = get(har, 'Extent');
    %     text(XX(1), Ynorm, Str, 'units', 'normalized', 'fontsize', 10, ...
    %         'fontweight', 'bold', 'horizontalalign', 'right', 'verticalalign', 'middle');
else,
    har = text(1, Ynorm, ['\leftarrow' Str], 'units', 'normalized', 'fontsize', 10, ...
        'fontweight', 'bold', 'horizontalalign', 'left', 'verticalalign', 'middle', varargin{:});
end

function local_cache(Y);
uix = Y.UniqueRecordingIndex;
GBdir = fullfile(processed_datadir,  '\JL\NNTP');
CFN = fullfile(GBdir, 'cycleCorrs');
load(CFN, 'CC');
ihit = find([CC.UniqueRecordingIndex]==uix);
CC(ihit) = Y;
save(CFN, 'CC', '-V6');











