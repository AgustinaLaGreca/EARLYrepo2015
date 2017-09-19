function P = JLpredictSpikes(Y, Thr, doPlot);
%  JLpredictSpikes - predict AP stats from JL binbeat data
%     JLpredictSpikes(Y, Thr, doPlot)
%
%     JLpredictSpikes(Y, nan) searches for Thr that reproduces # spikes

[Thr, doPlot] = arginDefaults('Thr/doPlot', nan, nargout<1);

Nspike = numel(Y.SPT); % true # spikes
if isnan(Thr), % find Thr
    local_match_Nspike('init', Y, Nspike);
    Thr = fminsearch(@local_match_Nspike, Y.SpontMean+sqrt(Y.SpontVar));
end

dt_bin = Y.bin_dt;
p_fire = SpikePredict(Thr, Y.bin_MeanWave, sqrt(Y.bin_VarWave));
AnDur = Y.NbeatCycle*Y.BeatCycle; % ms duration of analysis window
APrate = 1e3*p_fire/AnDur; % instantaneous rate spikes/s
pred_Nspike = sum(p_fire);
TotRate = 1e3*pred_Nspike/AnDur;
[dt_ipsi, Ipsi_rate] = fracLoopMean(dt_bin, APrate, Y.Period1);
[dt_contra, Contra_rate] = fracLoopMean(dt_bin, APrate, Y.Period2);
%dplot(dt_ipsi, Ipsi_rate);  xdplot(dt_contra, Contra_rate, 'r');  

timelag = Y.tstart-Y.preDur; % lag of analysis window re stimulus onset
pred_VSbin = local_VS(dt_bin, p_fire, Y.Fbeat, timelag);
pred_VSipsi = local_VS(dt_ipsi, Ipsi_rate, Y.Freq1, timelag);
pred_VScontra = local_VS(dt_contra, Contra_rate, Y.Freq2, timelag);
% data
Fbeat = Y.Fbeat; BeatCycle = 1e3/Fbeat;
Fipsi = Y.Freq1; Fcontra = Y.Freq2; 
IpsiCycle = Y.Period1; ContraCycle = Y.Period2;
Awin = Y.tstart + [0 AnDur]; % analysis window for cycle histograms
SPTstim = Y.SPTraw(betwixt(Y.SPTraw, Awin))-Y.preDur; % spike times during stim, re stim onset
Nspike = numel(SPTstim);
SPTbin = mod(SPTstim, BeatCycle); % SPTs wrapped over beat freq
SPTipsi = mod(SPTstim, IpsiCycle); % SPTs wrapped over ipsi freq
SPTcontra = mod(SPTstim, ContraCycle); % SPTs wrapped over contra freq
% ID
[ExpID, RecID, seriesID, iGerbil, icond, icell] ...
    = deal(Y.ExpID, Y.RecID, Y.seriesID, Y.iGerbil, Y.icond, Y.icell);
P = CollectInStruct(ExpID, RecID, seriesID, iGerbil, icond, icell, '-', ...
    Thr, dt_bin, p_fire, APrate, pred_Nspike, pred_VSbin, ...
    dt_ipsi, Ipsi_rate, pred_VSipsi, dt_contra, Contra_rate, pred_VScontra, '-', ...
    Fbeat, Fipsi, Fcontra, BeatCycle, IpsiCycle, ContraCycle, '-', ...
    Nspike, SPTstim ,SPTbin, SPTipsi, SPTcontra);

%=======================================
%============PLOT=======================
%=======================================
if ~doPlot, return; end
%NbinHisto = 15;
fh = figure;
set(fh,'units', 'normalized', 'position', [0.0219 0.655 0.977 0.247]);
figname(fh, [mfilename 'inst APrate']);
dplot(Y.bin_dt, APrate);
xlabel('Time (ms)');
ylabel('Inst firing rate (sp/s)');
title(sprintf('Total %0.0f spikes (%0.1f sp/s)', pred_Nspike, TotRate));
TracePlotInterface(fh)

fh2 = figure;
set(fh2,'units', 'normalized', 'position', [0.493 0.321 0.468 0.414])
figname(fh2, [mfilename 'Cycle Histo']);
% ipsi
h1=subplot(2,3,1);
Nbin = local_Chist(SPTstim, Y.Freq1);
title(sprintf('Ipsi %i Hz', round(Y.Freq1)));
ylabel('Measured # spikes')
h2=subplot(2,3,4);
binCenter = ((1:Nbin)-0.5)/Nbin;
BinHist = local_binnify(Ipsi_rate, Nbin);
BinHist = pred_Nspike*BinHist/sum(BinHist);
bar(binCenter, BinHist, 1);
local_uptext(sprintf('VS = %1.2f', abs(pred_VSipsi)));
ylabel('Predicted #spikes')
local_eqax(h1,h2)
% contra
h1=subplot(2,3,2);
Nbin = local_Chist(SPTstim, Y.Freq2);
title(sprintf('Contra %i Hz', round(Y.Freq2)));
h2=subplot(2,3,5);
binCenter = ((1:Nbin)-0.5)/Nbin;
BinHist = local_binnify(Contra_rate, Nbin);
BinHist = pred_Nspike*BinHist/sum(BinHist);
bar(binCenter, BinHist, 1);
ylim([min(ylim) mean(ylim)+0.6*diff(ylim)]);
local_uptext(sprintf('VS = %1.2f', abs(pred_VScontra)));
xlabel('Stim Phase (cycle)');
local_eqax(h1,h2)
% beat
h1=subplot(2,3,3);
Nbin = local_Chist(SPTstim, Y.Fbeat);
title(sprintf('Beat %i Hz', round(Y.Fbeat)));
h2=subplot(2,3,6);
binCenter = ((1:Nbin)-0.5)/Nbin;
BinHist = local_binnify(APrate, Nbin);
BinHist = pred_Nspike*BinHist/sum(BinHist);
bar(binCenter, BinHist, 1);
local_uptext(sprintf('VS = %1.2f', abs(pred_VSbin)));
local_eqax(h1,h2)
%dplot(dt, smoothen(APrate, 50, dt));



%================================================================
function Y = local_binnify(Y, Nbin);
N = numel(Y);
Nnew = Nbin*round(N/Nbin);
Y = interp1(linspace(0,1,N), Y(:), linspace(0,1,Nnew));
Y = reshape(Y, [Nnew/Nbin Nbin]);
Y = mean(Y,1).';

function VS = local_VS(dt, p_inst, Fana, tstart);
% dt in ms; p_inst arb scaling; Fana in Hz; tstart in ms
Cycle = 1e3/Fana; % analysis cycle in ms
Nsam = numel(p_inst);
NsamCycle = round(Cycle/dt);
NsamSel = NsamCycle*round(Nsam/NsamCycle); % select integer # cycles
if NsamSel<Nsam, NsamSel = NsamSel - NsamCycle; end
p_inst = p_inst(1:NsamSel);
Time = tstart + (0:NsamSel-1)*dt; % time re start of stimulus
Phase = Time/Cycle; % phase in cycles of the analysis freq
VS = sum(exp(2*pi*i*Phase(:)).*p_inst(:))/sum(p_inst(:));

function Mismatch = local_match_Nspike(Thr, Y, Nspike);
persistent YY NNspike
if isequal('init', Thr),
    YY = Y;
    NNspike = Nspike;
    return;
end
P = JLpredictSpikes(YY, Thr,0);
Mismatch = abs(P.pred_Nspike-NNspike);

function Nbin = local_Chist(Spt, Fana);
% Fana in Hz, Spt and Awin in ms
Nspike = numel(Spt);
if Nspike==0, return; end
R = abs(vectorstrength(Spt,Fana));
anaCycle = 1e3/Fana;
Spt = mod(Spt, anaCycle)/anaCycle;
Nbin = min(50,round(Nspike/10));
Nbin = max(Nbin, 10);
Ph = ((1:Nbin)-0.5)/Nbin;
hist(Spt, Ph);
local_uptext(sprintf('VS = %1.2f', R));

function local_uptext(Txt, varargin);
ylim([min(ylim) mean(ylim)+0.6*diff(ylim)]);
%text(mean(xlim)-0.3*diff(xlim), mean(ylim)+0.42*diff(ylim), Txt, varargin{:});
text(0.2, 0.93, Txt, 'units', 'normalized', varargin{:});

function local_eqax(h1,h2);
Y1 = get(h1,'ylim');
Y2 = get(h2,'ylim');
YL = [min(Y1(1),Y2(1)) max(Y1(2),Y2(2))];
set([h1 h2], 'ylim', YL);






