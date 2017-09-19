function [T, Tstats] = JLanova2(Uidx, doPlot);
% JLanova2 - two-ear anova of JL data
%   [T, Tstats] = JLanova2(S, doPlot);
%    
%   Tstats store in dbase D:\processed_data\JL\JLanova2\AnovaStats.dbase
%
%

[doPlot] = arginDefaults('doPlot', nargout<1);

S = JLdatastruct(Uidx); % specific info on stimulus etc
CFN = fullfile('D:\processed_data\JL\JLanova2', [mfilename '_' num2str(S.UniqueSeriesIndex)]);
DBFN = 'D:\processed_data\JL\JLanova2\AnovaStats.dbase';
init_dbase(DBFN, 'UniqueRecordingIndex', 'onlyifnew');
[T_Tstats CFN, CP] = getcache(CFN, Uidx);
if ~isempty(T_Tstats),
    [T, Tstats] = deal(T_Tstats.T, T_Tstats.Tstats);
else, % compute
    D = JLdbase(Uidx); % general info on the recording & stimulus
    W = JLwaveforms(Uidx); % waveforms and duration params
    T = JLcycleStats(Uidx); % cycle-based means & variances

    dt = W.dt; % ms sample period of detrended waveform
    freqI = S.Freq1; % ipsi freq
    freqC = S.Freq2; % contra freq

    Time = Xaxis(W.stimrec_detrended, dt, W.t_anaStart-W.t_stimOnset); % time of rec since stim onset
    NsamPerCycle = 40; % #bins in monaural cycle for 2-D plots
    iphaseI = floor(0.5+NsamPerCycle*rem(Time*1e-3*freqI, 1));
    iphaseC = floor(0.5+NsamPerCycle*rem(Time*1e-3*freqC, 1));
    iphaseIC = iphaseI + NsamPerCycle*iphaseC;
    % hist(NsamPerCycle*iphase1+iphase2,0:NsamPerCycle^2+1)
    %[P,T,STATS,TERMS] = anovan(rec,{iphase1,iphase2});

    % ====spike times====
    SPT = W.SPTraw-W.t_stimOnset; % spike times expressed re stim onset
    SPT = SPT(betwixt(SPT, min(Time), max(Time))); % only consider spikes occurring within analysis window
    NspikeStim = numel(SPT);
    Spike_phase_I = rem(SPT*1e-3*freqI, 1);
    Spike_phase_C = rem(SPT*1e-3*freqC, 1);
    isp_phaseI = floor(NsamPerCycle*Spike_phase_I);
    isp_phaseC = floor(NsamPerCycle*Spike_phase_C);

    % visit each possible combination of monaural (binned) phases and compute
    % mean and variance over that bin, and the count the # spikes falling
    % within that bin.
    [mean2d, var2d, Nspike2d] = deal(zeros(NsamPerCycle));
    for i1=0:NsamPerCycle-1,
        for i2=0:NsamPerCycle-1,
            %qcell = (iphaseI==i1)&(iphaseC==i2);
            qcell = iphaseIC == i1+NsamPerCycle*i2;
            mean2d(i1+1,i2+1) = nanmean(W.stimrec_detrended(qcell));
            var2d(i1+1,i2+1) = nanvar(W.stimrec_detrended(qcell));
            Nspike2d(i1+1,i2+1) = sum((isp_phaseI==i1)&(isp_phaseC==i2));
        end
    end
    phi_mon = (0:NsamPerCycle-1)/NsamPerCycle;
    [Phase_C, Phase_I]= meshgrid(phi_mon,phi_mon);
    Phase_B = deciRound(mod(Phase_C - Phase_I,1),4); % interaural phase, rounded to eliminate rounding errors :)
    phi_bin = unique(Phase_B(:)); % unique values of interaural phases
    for ii=1:numel(phi_bin), % mean of values having interaural phase UDphi(ii)
        binmean(ii) = mean(mean2d(Phase_B==phi_bin(ii)));
    end
    % attempt to remove incidental nans from mean2d
    mean2d = local_removeNans(mean2d);
    var2d = local_removeNans(var2d);

    ipsimean = mean(mean2d,2);
    contramean = mean(mean2d,1);
    % evaluate lag between subthr input & spikes.
    if numel(SPT)>0,
        XX = real(ifft2(fft2(Nspike2d).*conj(fft2(mean2d)))); % circular 2D corr fnc XX(1,1) is zero/zero lag
        [dum, imax] = max(XX(:));
        lag_ipsi = Xaxis(1:NsamPerCycle, 1e3/freqI/NsamPerCycle);
        lag_contra = Xaxis(1:NsamPerCycle, 1e3/freqC/NsamPerCycle);
        [lag_contra, lag_ipsi] = meshgrid(lag_contra, lag_ipsi);
        [APlagI, APlagC] = deal(lag_ipsi(imax), lag_contra(imax));
    else,
        [APlagI, APlagC] = deal(nan);
    end
    % for each spike, determine the mean & var of the subthr potential at its
    % moment in the stimulus, while taking into account the lags just computed.
    SPT_i = SPT - APlagI; % spike times corrected for ipsi lag
    SPT_c = SPT - APlagC;  %ditto contra
    phi_i = rem(SPT_i*1e-3*freqI, 1); % Ipsi phase at the corrected spike times
    phi_c = rem(SPT_c*1e-3*freqC, 1); % ditto contra
    isp_phaseI = 1+floor(NsamPerCycle*phi_i);
    isp_phaseC = 1+floor(NsamPerCycle*phi_c);
    mean_at_spt = [];
    var_at_spt = [];
    for ii=1:numel(SPT), % mean & var potential at the corrected spike times
        mean_at_spt(ii) = mean2d(isp_phaseI(ii), isp_phaseC(ii));
        var_at_spt(ii) = var2d(isp_phaseI(ii), isp_phaseC(ii));
    end
    % from this, compute spike prob
    Vinst = linspace(min(mean2d(:)), max(mean2d(:)), 15);
    histV = hist(mean2d(:), Vinst);
    DurPerCell = 1e-3*T.AnaDur/NsamPerCycle.^2; % s duration represented by one cell in 2D grid
    hist_mean_at_spt = hist(mean_at_spt, Vinst);
    hist_var_at_spt = hist(var_at_spt, Vinst);
    cond_spike_prob = hist_mean_at_spt./(DurPerCell*histV);

    linpred = SameSize(ipsimean-0.5*mean(ipsimean), contramean) + SameSize(contramean-0.5*mean(contramean), ipsimean);
    UniqueRecordingIndex = Uidx;
    VarAcc = VarAccount(mean2d,linpred);
    SPTreStimOnset = SPT;
    T = CollectInStruct(UniqueRecordingIndex, '-anova2', NsamPerCycle, ...
        '-mean_and_var2D', mean2d, var2d, phi_mon, ipsimean, contramean, phi_bin, binmean, linpred, VarAcc, ...
        '-spike2D', SPTreStimOnset, Nspike2d, Spike_phase_I, Spike_phase_C ,APlagI, APlagC, ...
        '-spikeprob', NspikeStim, mean_at_spt, var_at_spt, Vinst, histV, hist_mean_at_spt, hist_var_at_spt, cond_spike_prob);
    % Single-number struct
    Anova_VarAcc = VarAcc;
    Tstats = CollectInStruct(UniqueRecordingIndex, '-Anova2_stats', Anova_VarAcc, APlagI, APlagC);
    % R = collectInStruct(dt, rec);
    % R = structjoin(R, '-', S);
    % caching
    putcache(CFN, 100, CP, CollectInStruct(T, Tstats));
    % primitive caching derived from JL_NTTP
    local_cache(T, Tstats);
end
Add2dbase(DBFN, Tstats);

%S = JLdatastruct(Uidx); % standardized recording info

if doPlot,
    local_plot(T,S);
end

%=========================================================
%=========================================================
%=========================================================

function ipeak = local_findpeak(X);
% peak of circular waveform
X = (X(:)-min(X(:))).^4;
Nsam = numel(X);
Ph = (0:Nsam-1).'/Nsam;
phaseFac = exp(2*pi*i*Ph);
PeakPhase = angle(sum(phaseFac.*X))/2/pi;
PeakPhase = mod(PeakPhase,1);
ipeak = 1+floor(PeakPhase*Nsam);

function  [spikelag_ipsi, spikelag_contra] = local_lag(Nspike, dt, ipsimean, contramean, S);
if sum(Nspike(:))==0,
    [spikelag_ipsi, spikelag_contra] = deal(nan);
    return;
end
Nsam = numel(ipsimean);
nsp1= sum(Nspike,2);
nsp2= sum(Nspike,1).';
plot(nsp1)
xplot(nsp2,'r')
xplot(ipsimean,'b:')
xplot(10*ipsimean,'b:')
xplot(10*contramean,'r:')


function local_plot(T,S);
if isempty(S), S = JLdatastruct(T.UniqueRecordingIndex); end
Phase = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
imaxI = local_findpeak(T.ipsimean);
imaxC = local_findpeak(T.contramean);
nshift = round(T.NsamPerCycle/2)-[imaxI,imaxC];
nshift = [0 0];%*nshift;

mean2d = circshift(T.mean2d, nshift);
linpred = circshift(T.linpred, nshift);
var2d = circshift(T.var2d, nshift);
Nspike = circshift(T.Nspike2d, nshift);
Spike_phase_C = rem(T.Spike_phase_C+2+nshift(1)/T.NsamPerCycle, 1);
Spike_phase_I = rem(T.Spike_phase_I+2+nshift(2)/T.NsamPerCycle, 1);

set(figure,'units', 'normalized', 'position', [0 0.183 0.864 0.722], 'PaperOrientation', 'landscape')
%========mean waves vs phase===========
% DT_I = 1e3/S.Freq1/T.NsamPerCycle;
% DT_C = 1e3/S.Freq1/T.NsamPerCycle;
axes('position', [0.2183    0.1767    0.1172    0.2121]);
plot(T.phi_mon,T.ipsimean(:), 'b');
xplot(T.phi_mon,T.contramean(:), 'r');
xplot(T.phi_bin,T.binmean(:), 'k');
xlabel('Phase (cycle)');
ylabel('Avg Potential (mV)');
legend({'I' 'C' 'B'}, 'location', 'NorthOutside')
set(gca, 'position', [0.2328    0.1496    0.1172    0.2121])
%==========surface & contour plots of mean & variance============
VV = [ 0.5592    0.8290         0   -0.6941
    -0.6351    0.4284    0.6428   -0.2180
    -0.5329    0.3594   -0.7660    9.1300
    0         0         0    1.0000];

subplot(4,3,[8 11]);
surface(Phase, Phase, mean2d);
view(VV);
zlabel('Mean Potential (mV)')
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')

subplot(4,3,[2 5]);
surface(Phase, Phase, var2d);
view(VV);
zlabel('Variance (mV^2)')
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')

subplot(4,3,[9 12]);
contourf(Phase, Phase, mean2d);
axis square
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')
xplot(Spike_phase_C, Spike_phase_I, '.w');
title(['spike lags [' num2str(round(T.APlagI*1e3)) ', '  num2str(round(T.APlagC*1e3)) '] us'])

subplot(4,3,[3 6]);
contourf(Phase, Phase, var2d);
axis square
irm = strfind(S.TTT,'4 Hz');
ht = title(sprintf('exp %d  cell %d:  %s, %d Hz, %d dB SPL (series %d).', ...
    S.iexp, S.icell, S.chan, S.freq, S.SPL, S.iseries), 'fontsize' ,14);
tpos = get(ht,'position');
set(ht,'position', [tpos(1)-0.5*diff(tpos([1 3])) tpos(2)]);
% set(gcf,'paperpositionmode', 'auto')
%xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')


%=======deviation from I+C prediction =======
Deviation = T.mean2d - T.linpred;
%Deviation = Deviation./sqrt(mean(var2d(:)));
axes('position', [0.0859    0.6760    0.1292    0.2872]);
contourf(Phase, Phase, Deviation);
axis square;
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')
title('Deviation from linearity');
hc = colorbar;
set(hc,'yaxisloc', 'right');
ylabel(hc, 'Deviation (mV)');
set(gca, 'position', [0.0859    0.6760    0.1292    0.2872]);
% VarAcc
axes('position', [0.0522    0.4061    0.1731    0.2583]);
plot(T.linpred, T.mean2d, '.');
xplot(xlim, xlim, 'k', 'linewidth', 2);
text(0.1,0.8, ['VarAcc = ' num2str(deciRound(T.VarAcc,4)) ' %'], 'units', 'normalized');
xlabel('V_{pred} (mV)');
ylabel('V_{rec} (mV)');
%set(gcf,'paperpositionmode', 'auto', 'renderer', 'zbuffer')

%====firing prob against inst potential, corrected for lag============
axes('units', 'normaliz', 'position', [0.0530    0.1046    0.1300    0.2180]);
plot(T.Vinst, T.cond_spike_prob, '*-');
xlabel('Potential (mV)');
ylabel('Instant. firing rate (1/s)');
%set(gcf,'paperpositionmode', 'auto', 'renderer', 'zbuffer')
text(0.1, 0.85, [num2str(T.NspikeStim) ' APs'], 'units', 'normalized')

%=============================================================
function X = local_removeNans(X);
% replace nans by average btw neighbors
inan = find(isnan(X(:)));
iplus = 1+rem(inan,numel(X));
imin = 1+rem(inan-2,numel(X));
X(inan) = (X(imin)+X(iplus))/2;

function local_cache(A, astats);
% primitive caching derived from JL_NTTP
uix = A.UniqueRecordingIndex;
GBdir = fullfile(processed_datadir,  '\JL\NNTP');
load(fullfile(GBdir, 'Anova2_Stats'), 'Astats');
ihit = find([Astats.UniqueRecordingIndex]==uix);
Astats(ihit) = astats;
save(fullfile(GBdir, 'Anova2_Stats'), 'Astats', '-V6');
save(fullfile(GBdir, ['Anova2_details_' num2str(uix)]), 'A');




