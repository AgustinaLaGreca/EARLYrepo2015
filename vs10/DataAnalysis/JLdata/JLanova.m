function [T, R] = JLanova(S, varargin);
% JLanova - two-ear anova of JL data
%   T = JLanova(S);
%
%   [T, R] = JLanova(S); 
%   R is resampled recording
%
%   JLanova('plot', T)


if isequal('plot', S),
    local_plot(varargin{1}, []);
    return;
end


S = JLdatastruct(S); % standardized recording info
CFN = fullfile(processed_datadir, 'JL', 'JLanova', [S.ExpID '__' S.RecID '.cache']);
[T, CFN, Param] = getcache(CFN, S.UniqueRecordingIndex);
%T = [] % switch caching off

if isempty(T) || nargout>1, % compute T
    preDur = 500; % ms pre-stim by convention
    postDur = 500; % ms post stim after which to analyze spont
    tstartSpont = preDur+S.burstDur+postDur; % evaluate spont activity starting 500 ms after stim offset
    Ttrans = S.Ttrans; % transition time after onset, excluded from analysis

    Uidx = S.UniqueRecordingIndex;
    %S = JLbeatList(S.ExpID, S.icell);
    S = S([S.UniqueRecordingIndex]==Uidx);

    [D, DS, L]=readTKABF(S);
    dt = D.dt_ms;
    rec = D.AD(1).samples; % entire recording, including pre- and post-stim parts
    %rec = D.AD(3).samples + D.AD(4).samples;  %DEBUG  stimuli of the two ears summed
    NsamTot = D.NsamPerSweep/D.Nchan;
    cutWin = [S.APwindow_start S.APwindow_end];
    [rec, SPTraw, Tsnip, Snip, APslope] = ...
        APtruncate2(rec, [S.APthrSlope, S.CutoutThrSlope], dt, cutWin, 1); % last 1: replace APs by nans
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

    % resample to make beatcycle contain integer # samples
    oldNsam = numel(rec);
    Nsam = NbeatCycle*round(oldNsam/NbeatCycle);
    rec = interp1(linspace(0,1,oldNsam), rec, linspace(0,1,Nsam));
    dt_old = dt;
    dt = dt_old*oldNsam/Nsam;

    Time = Xaxis(rec, dt, Ttrans); % time since stim onset
    totDur = numel(rec)*dt; % duration of rec
    Period1 = 1e3/S.Freq1;
    NsamPerCycle = round(Period1/dt);
    NsamPerCycle = 40; % #bins in monaural cycle for 2-D plots
    iphaseI = floor(NsamPerCycle*rem(Time*1e-3*S.Freq1, 1));
    iphaseC = floor(NsamPerCycle*rem(Time*1e-3*S.Freq2, 1));
    % hist(NsamPerCycle*iphase1+iphase2,0:NsamPerCycle^2+1)
    %[P,T,STATS,TERMS] = anovan(rec,{iphase1,iphase2});

    % ====spike times====
    SPT = SPTraw-preDur; % spike times expressed re stim onset
    SPT = SPT(betwixt(SPT, min(Time), max(Time))); % only consider spikes occurring within analysis window
    Spike_phase_I = rem(SPT*1e-3*S.Freq1, 1);
    Spike_phase_C = rem(SPT*1e-3*S.Freq2, 1);
    isp_phaseI = floor(NsamPerCycle*Spike_phase_I);
    isp_phaseC = floor(NsamPerCycle*Spike_phase_C);

    % visit each possible combination of monaural (binned) phases and compute
    % mean and variance over that bin, and the count the # spikes falling
    % within that bin.
    [mean2d, var2d] = deal(zeros(NsamPerCycle));
    for i1=0:NsamPerCycle-1,
        for i2=0:NsamPerCycle-1,
            qcell = (iphaseI==i1)&(iphaseC==i2);
            mean2d(i1+1,i2+1) = nanmean(rec(qcell));
            var2d(i1+1,i2+1) = nanvar(rec(qcell));
            Nspike(i1+1,i2+1) = sum((isp_phaseI==i1)&(isp_phaseC==i2));
        end
    end
    
    ipsimean = mean(mean2d,2);
    contramean = mean(mean2d,1);
    %[spikelag_ipsi, spikelag_contra] = local_lag(Nspike, dt, ipsimean, contramean, S);
    linpred = SameSize(ipsimean-0.5*mean(ipsimean), contramean) + SameSize(contramean-0.5*mean(contramean), ipsimean);
    UniqueRecordingIndex = S.UniqueRecordingIndex;
    VarAcc = VarAccount(mean2d,linpred),4;
    T = CollectInStruct(UniqueRecordingIndex, NsamPerCycle, dt, '-', ...
        mean2d, var2d, ipsimean, contramean, linpred, VarAcc ,'-', ...
        Nspike, Spike_phase_I, Spike_phase_C);
    putcache(CFN, 40, Param, T); 
    R = CollectInStruct(dt, rec);
    R = structJoin(R, '-', S);
else,
    S = [];
end % if ~isempty(T)


local_plot(T,S);

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


function local_plot(T, S);
if isempty(S), S = JLdatastruct(T.UniqueRecordingIndex); end
Phase = linspace(-0.5,0.5,T.NsamPerCycle);
imaxI = local_findpeak(T.ipsimean);
imaxC = local_findpeak(T.contramean);
nshift = round(T.NsamPerCycle/2)-[imaxI,imaxC];
nshift = [0 0];%*nshift;

mean2d = circshift(T.mean2d, nshift);
linpred = circshift(T.linpred, nshift);
var2d = circshift(T.var2d, nshift);
Nspike =  circshift(T.Nspike, nshift);
Spike_phase_C = rem(T.Spike_phase_C+2+nshift(1)/T.NsamPerCycle, 1)-0.5;
Spike_phase_I = rem(T.Spike_phase_I+2+nshift(2)/T.NsamPerCycle, 1)-0.5;

%========mean wave vs time===========
set(figure,'units', 'normalized', 'position', [0.0328 0.667 0.255 0.224]);
plot(Phase(:),T.ipsimean(:), 'b');
xplot(Phase(:),T.contramean(:), 'r');

%==========surface & contour plots of mean & variance============
set(figure,'units', 'normalized', 'position', [0.254 0.173 0.695 0.751])
VV = [ 0.5592    0.8290         0   -0.6941
    -0.6351    0.4284    0.6428   -0.2180
    -0.5329    0.3594   -0.7660    9.1300
    0         0         0    1.0000];

subplot(2,2,3);
surface(Phase, Phase, mean2d);
view(VV);
zlabel('Mean Potential (mV)')
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')

subplot(2,2,1);
surface(Phase, Phase, var2d);
view(VV);
zlabel('Variance (mV^2)')
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')

subplot(2,2,4);
contourf(Phase, Phase, mean2d);
axis square
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')
xplot(Spike_phase_C, Spike_phase_I, '.w');

subplot(2,2,2);
contourf(Phase, Phase, var2d);
axis square
irm = strfind(S.TTT,'4 Hz');
ht = title(S.TTT(1:irm-1), 'fontsize' ,14);
tpos = get(ht,'position');
set(ht,'position', [tpos(1)-0.5*diff(tpos([1 3])) tpos(2)]);
set(gcf,'paperpositionmode', 'auto')
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')


%=======deviation from I+C prediction =======
set(figure,'units', 'normalized', 'position', [0.0328 0.102 0.251 0.525])
Deviation = T.mean2d - T.linpred;
%Deviation = Deviation./sqrt(mean(var2d(:)));
subplot(2,1,1);
contourf(Phase, Phase, Deviation);
axis square;
xlabel('Contra phase (cycle)')
ylabel('Ipsi phase (cycle)')
title('Deviation from linearity');
hc = colorbar;
set(hc,'yaxisloc', 'right');
ylabel(hc, 'Deviation (mV)');
subplot(2,1,2);
plot(T.linpred, T.mean2d, '.');
xplot(xlim, xlim, 'k', 'linewidth', 2);
text(0.1,0.8, ['VarAcc = ' num2str(deciRound(T.VarAcc,4)) ' %'], 'units', 'normalized');
xlabel('V_{pred} (mV)');
ylabel('V_{rec} (mV)');
set(gcf,'paperpositionmode', 'auto', 'renderer', 'zbuffer')
