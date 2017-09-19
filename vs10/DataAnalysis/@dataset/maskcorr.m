function [AllCorr, CondMetrics] = maskcorr(D,D2);
% dataset/maskcorr - autocorrelogram analysis of MASK (tone+noise) data
%    S = maskcorr(DS) produces across-condition correlograms of MASK dataset DS
%  
%    maskcorr(DS1, DS2) produces across-dataset correlograms.
%
%    See also cyclehisto, SPTcorr.

MaxLag = 7; % ms max lag
BinWidth = nan; % ms bin width; depends on tone freq (see below)
minISI = 0.5; % minimum inter-spike-interval
rateAnDur = 400; % ms before & after tone onset for rate analysis 
if ~isequal('MASK', upper(stimtype(D))),
    error('First dataset is not a MASK dataset.');
end
ExpID = name(D.ID.Experiment);
iRec = D.ID.iDataset;

if isnan(BinWidth),
    if D.Stim.ToneFreq<=1000,
        BinWidth = 1e3/D.Stim.ToneFreq/10; % 1/10 of a cycle
    else,
        BinWidth = 0.1;
    end
end


set(figure,'units', 'normalized', 'position', [0 0.0267 1 0.91]);
fh = gcf;
[MinY, MaxY] = deal(0);
if nargin==1,
    SingleDS = 1;
    SPT = spiketimes(D,1,'', minISI);
    Ncond = size(SPT,1);
    [NoiseWin ToneWin] = local_Awin(D);
    for icond = 1:Ncond,
        spt_i = SPT(icond, :);
        warp1 = D.Stim.Warp(icond);
        for jcond = icond:Ncond,
            warp2 = D.Stim.Warp(jcond);
            spt_j = SPT(jcond, :);
            [SacN Tau SacVarN] = local_corr(icond==jcond, AnWin(spt_i, NoiseWin), AnWin(spt_j, NoiseWin), MaxLag, BinWidth, diff(NoiseWin));
            [SacT Tau SacVarT] = local_corr(icond==jcond, AnWin(spt_i, ToneWin), AnWin(spt_j, ToneWin), MaxLag, BinWidth, diff(ToneWin));
            ha(icond, jcond) = subplot(Ncond, Ncond, jcond+Ncond*(Ncond-icond), 'visible', 'off');
            totDiff = mean((sqrt(SacN)-sqrt(SacT)).^2); % "sum of squares" of sqrt(SacX) differences
            [dum, ipeak1] = min(abs(Tau));
            [dum, ipeak2] = min(abs(Tau-1e3/D.Stim.ToneFreq));
            peak1Diff = sqrt(SacT(ipeak1))-sqrt(SacN(ipeak1)); 
            peak2Diff = sqrt(SacT(ipeak2))-sqrt(SacN(ipeak2)); 
            plot(Tau, sqrt(SacN), 'linewidth',2, 'visible', 'off');
            xplot(Tau(ipeak1), sqrt(SacN(ipeak1)), 'p', 'visible', 'off');
            xplot(Tau(ipeak2), sqrt(SacN(ipeak2)), 'o', 'visible', 'off');
            DsqN = 0.5*sqrt(SacVarN./SacN); % std of sqrt(SacN)
            local_stdplot(Tau,sqrt(SacN),DsqN, [0 0 1], 2);
            xplot(Tau, sqrt(SacT), 'r', 'linewidth',2, 'visible', 'off');
            SacDiff = sqrt(SacT)-sqrt(SacN);
            fitcos = cos(2*pi*Tau(:)*1e-3*D.Stim.ToneFreq);
            CosfitSacDiff = corr(SacDiff(:),fitcos);
            xplot(Tau, SacDiff, 'color', [0 0.75 0], 'visible', 'off');
            xplot(Tau, 0*SacDiff, ':', 'color', [0 0.75 0], 'visible', 'off');
            xplot([min(Tau) max(Tau)], [1 1], 'k', 'visible', 'off');
            Lbl = trimspace(num2str([icond jcond]));
            Lbl = strrep(Lbl ,' ', '/');
            text(0.1, 0.85, Lbl, 'units', 'normalized', 'visible', 'off');
            pkstrN = num2str(deciRound(max(SacN),3));
            pkstrT = num2str(deciRound(max(SacT),3));
            text(0.7, 0.85, pkstrN, 'units', 'normalized', 'color', 'b', 'visible', 'off');
            text(0.7, 0.7, pkstrT, 'units', 'normalized', 'color', 'r', 'visible', 'off');
            set(gca,'xtick', [], 'ytick', []);
            set(gca,'position', get(gca, 'position').*[1 1 1.2 1.3]);
            MaxY = max(MaxY, max([max(sqrt(SacN)+DsqN) max(sqrt(SacT))]));
            MinY = min(MinY, min(SacDiff));
            set(gca, 'userdata', [icond, jcond]); % used by callback fcns like local_fano
            AllCorr(icond,jcond) = CollectInStruct(ExpID, iRec, icond, jcond, warp1, warp2,'-', ...
                Tau, SacN, SacT, SacVarN, SacVarT, '-', ipeak1, ipeak2, '-', ...
                totDiff, peak1Diff, peak2Diff, CosfitSacDiff);
            AllCorr(jcond,icond) = AllCorr(icond,jcond);
            [AllCorr(jcond,icond).icond, AllCorr(jcond,icond).jcond] = swap(AllCorr(jcond,icond).icond, AllCorr(jcond,icond).jcond);
            [AllCorr(jcond,icond).warp1, AllCorr(jcond,icond).warp2] = swap(AllCorr(jcond,icond).warp1, AllCorr(jcond,icond).warp2);
        end % for jcond
    end % for icond
    setGUIdata(fh, 'AllCorr', AllCorr);
else, % across datasets
    if ~isequal('MASK', stimtype(D)),
        error('Second dataset is not a MASK dataset.');
    end
    error NYI
end

% linkaxes(ha(ha>0), 'xy'); this severly bogs down closing the figure 
%linkaxes(ha(ha>0), 'off');
for ii=1:Ncond,
    axes(ha(1,ii));
    set(gca,'xtick', -20:5:20);
    xlabel(num2str(D.Stim.Warp(ii)), 'fontsize', 14, 'fontweight', 'bold');
    axes(ha(ii,end));
    set(gca,'ytick', -20:1:20, 'YAxisLocation', 'right');
    ylabel(num2str(D.Stim.Warp(ii)), 'fontsize', 14, 'fontweight', 'bold');
end
rha = ha(ha>0);
for ii=1:numel(rha),
    xlim(gca, MaxLag*[-1 1]);
    axes(rha(ii));
    ylim([MinY MaxY]); % manual "link axes"; much faster when closing window
end
MinY, MaxY
set(findobj(fh, 'type', 'line'), 'visible', 'on');
set(findobj(fh, 'type', 'axes'), 'visible', 'on');
set(findobj(fh, 'type', 'patch'), 'visible', 'on');
% figure settings for printing
set(gcf,'PaperOrientation', 'landscape', 'paperunits', 'inches', 'paperpos',[-0.4 -0.3 11.5 9]);

%=================================
%==========extra plots============
%=================================
% rate plot
figure(fh);
hr  = axes('units', 'normalized', 'position', [0.1 0.8 0.15 0.12]);
CondMetrics = local_rateplot(hr, D, minISI, rateAnDur);
% vector strength
hv  = axes('units', 'normalized', 'position', [0.1 0.65 0.15 0.12]);
VS = local_vsplot(hv, D, minISI, rateAnDur ,AllCorr);
CondMetrics = structJoin(CondMetrics, VS);
% SSQ diff
figure(fh);
hs  = axes('units', 'normalized', 'position', [0.3 0.8 0.15 0.12]);
local_ssqplot(hs, AllCorr);
% Peak
figure(fh);
hp  = axes('units', 'normalized', 'position', [0.3 0.65 0.15 0.12]);
local_peakplot(hp, AllCorr);
% CF from SAC spec
hc = axes('units', 'normalized', 'position', [0.1 0.5 0.15 0.12]);
local_SACspec(D, 0.15);

% display info about dataset
ht = axes('units', 'normalized', 'position', [0.1 0.35 0.25 .15]);
set(ht, 'visib' ,'off');
if nargin ==1,
    St = D.Stim;
    Pres = St.Presentation;
    Str = strvcat([name(D.ID.Experiment) ' rec ' num2str(D.ID.iDataset) '   <' IDstring(D) '>']);
    DurStr = [num2str(Pres.Ncond) ' x ' num2str(Pres.Nrep) ' x ' num2str(St.ISI) ' ms'];
    SPLstr = ['Noise: ' num2str(St.MnoiseSPL) ' ' St.MnoiseSPLUnit];
    ToneStr = ['Tone: ' num2str(St.ToneFreq) ' Hz  S/N = ' num2str(St.ToneSPL)  ' dB'];
    [xmin, xmax] = minmax(Pres.X.PlotVal(:)); [xmin, xmax] = deal(deciRound(xmin), deciRound(xmax));
    %ParStr = [Pres.X.ParName ' = ' num2str(xmin) ':' num2str(xmax) ' ' Pres.X.ParUnit];
    Str = strvcat(Str, SPLstr, ToneStr, DurStr);
else,
end
text(0.5, 0.5, Str, 'fontsize', 11, 'horizontalalign', 'center', 'verticalalign', 'middle');

% store axes sizes
ha = findobj(fh, 'type', 'axes');
for ii=1:numel(ha),
    AxProp(ii).handle = ha(ii);
    AxProp(ii).position = get(ha(ii), 'position');
end
setGUIdata(fh, 'AxProp', AxProp);
set(fh,'keypressfcn', @local_keypress);

%==========================
function [Sac Tau SacVar] = local_corr(Nodiag, SPT1, SPT2, MaxLag, BinWidth, Dur);
% cell arrays SPT1 & SPT2. cells are reps.
if Nodiag, % within dataset and within condition; exclude diag
    [SacVar Tau Sac] = SPTcorrvar(SPT1, 'nodiag', MaxLag, BinWidth, Dur, 'DriesNorm');
else, % across conditions and/or datasets
    [SacVar Tau Sac] = SPTcorrvar(SPT1, SPT2, MaxLag, BinWidth, Dur, 'DriesNorm');
end


function [NoiseWin ToneWin] = local_Awin(ds);
if ds.Stim.ToneOnset>0, % noise leads tone
    NoiseWin = [0 ds.Stim.ToneOnset];
    ToneWin = ds.Stim.ToneOnset+[0 ds.Stim.ToneDur];
else, % tone leads noise
    NoiseWin = -ds.Stim.ToneOnset+ + [0 ds.Stim.MnoiseBurstDur];
    ToneWin = [0 ds.Stim.ToneDur];
end

function Rate = local_rateplot(ha, D, minISI, AnDur);
SPT = spiketimes(D, 1, 'no-unwarp', minISI);
TimeFactor = 2.^(-D.Stim.Warp);
Ncond = numel(TimeFactor);
ToneOnset = D.Stim.ToneOnset;
Nrep = NrepRec(D);
Nrep = min(Nrep(:));
for icond = 1:Ncond,
    tfac = TimeFactor(icond);
    sptN = AnWin(SPT(icond,1:Nrep), tfac*(ToneOnset+[-AnDur 0]));
    sptT = AnWin(SPT(icond,:), tfac*(ToneOnset+[0 AnDur]));
    RateN(icond,:) = 1e3*cellfun(@numel,sptN)/(tfac*AnDur);
    RateT(icond,:) = 1e3*cellfun(@numel,sptT)/(tfac*AnDur);
end
X = D.Stim.Presentation.X;
DX = mean(diff(X.PlotVal));
warp = X.PlotVal;
% mean rates and std
meanRateN = mean(RateN,2);
stdRateN = std(RateN,[], 2);
meanRateT = mean(RateT,2);
stdRateT = std(RateT,[], 2);
%plot
axes(ha);
hN = errorbar(-0.1*DX + warp, meanRateN, stdRateN);
hold on
hT = errorbar(0.1*DX + warp, meanRateT, stdRateT);
set(hT,'color', 'r');
xplot(-0.1*DX + warp, meanRateN, 'b-o');
xplot(0.1*DX + warp, meanRateT, 'r-s');
ylim([0 max(ylim)]);
xlim([min(X.PlotVal)-0.5*DX max(warp)+0.5*DX])
%xlabel('Warp (Oct)');
ylabel('Firing rate sp/s');
legend({'Noise' 'Noise+Tone'}, 'fontsize', 6, 'location', 'best');
legend(ha, 'boxoff');
set(ha,'xticklabel',{});
Rate = CollectInStruct(warp, meanRateN, stdRateN, meanRateT, stdRateT);

function Y = local_vsplot(ha, D, minISI, AnDur, AC);
SPT = spiketimes(D, 1, 'no-unwarp', minISI);
TimeFactor = 2.^(-D.Stim.Warp);
Ncond = numel(TimeFactor);
ToneOnset = D.Stim.ToneOnset;
Nrep = NrepRec(D);
Nrep = min(Nrep(:));
Freq = D.Stim.ToneFreq;
for icond = 1:Ncond,
    tfac = TimeFactor(icond);
    spt = AnWin(SPT(icond,1:Nrep), tfac*(ToneOnset+[0 AnDur])); % select tone part
    VS(:,icond)  = vectorstrength(spt,Freq/tfac);
    [dum, Alpha(icond,1)] = vectorstrength([spt{:}],Freq/tfac); % evaluate Rayleigh significance from pooled reps
end
VS = abs(VS).'; % columns are now reps, rows, conditions
X = D.Stim.Presentation.X;
DX = mean(diff(X.PlotVal));
meanVS = mean(VS,2);
stdVS = std(VS,[], 2);
axes(ha);
hN = errorbar(-0.1*DX + X.PlotVal, meanVS, stdVS);
set(hN,'color','r');
xplot(-0.1*DX + X.PlotVal, meanVS, 'r-o', 'markerfacecolor', 'w');
xplot(-0.1*DX + X.PlotVal, pmask(Alpha<=0.001)+meanVS, 'r-o', 'markerfacecolor', 'r');
ylim([0 max(ylim)]);
xlim([min(X.PlotVal)-0.5*DX max(X.PlotVal)+0.5*DX])
%xlabel('Warp (Oct)');
ylabel('Vector strength');
ylim([0 1]);
% cosfit
ACd = AC([AC.icond]==[AC.jcond]);
xplot([ACd.warp1], [ACd.CosfitSacDiff], '-d', 'color', 0.75*[0 1 0]);
Y = CollectInStruct(meanVS, stdVS, Alpha);

function local_zoom(ha, zfac);
pos = get(ha,'position');
pos(3:4) = pos(3:4)*zfac;
if sum(pos([1 3]))>1, pos(1)=0.95-pos(3); end
if sum(pos([2 4]))>1, pos(2)=0.95-pos(4); end
set(ha,'position',pos);
% move to foreground
figh = parentfigh(ha);
ch = get(figh,'children');
ihit = find(ch==ha);
ch = ch([ihit 1:ihit-1  ihit+1:end]);
set(figh, 'children', ch);
drawnow;

function local_restoresize(figh, ah);
AxProp = getGUIdata(figh, 'AxProp');
if nargin>1, % select axes
    AxProp = AxProp([AxProp.handle]==ah);
end
for ii=1:numel(AxProp),
    set(AxProp(ii).handle, 'position', AxProp(ii).position);
end

function local_fano(Src, ah);
uc = get(ah, 'userdata');
if isempty(uc), return; end
[icond, jcond] = DealElements(uc);
AllCorr = getGUIdata(parentfigh(Src), 'AllCorr');
CC = AllCorr(icond, jcond);
figure;
plot(CC.SacN, CC.SacVarN, '.');
xplot(CC.SacT, CC.SacVarT, 'r.');
xlim([0 max(xlim)]); ylim([0 max(ylim)]);
xplot(xlim, xlim, 'k'); xplot(ylim, ylim, 'k');
xlabel('Corr mean');
ylabel('Corr var');


function local_keypress(Src, Ev);
%Ev
switch Ev.Key,
    case 'space';
        local_restoresize(Src);
    case 'escape';
        commandwindow;
    case 'delete';
        close(parentfigh(Src));
    case 'equal'; % +
        local_zoom(gca,sqrt(2));
    case 'hyphen'; % -
        local_zoom(gca,sqrt(1/2));
    case '0'; % -
        local_restoresize(Src,gca);
    case 'v'; % -
        local_fano(Src,gca);
    case {'h', '/'}; % -
        hStr = strvcat('space: restore all graph sizes', ...
            'escape: to MATLAB command window', ...
            'delete: close figure window', ...
            '     +: zoom in current graph', ...
            '     -: zoom out current graph', ...
            '     0: restore size of current graph', ...
            '     v: var-versus-mean analysis', ...
            '     h: help');
        figure; 
        axes('visible', 'off')
        text(0.5, 0.5, hStr, 'horizontalalign', 'center', 'verticalalign', 'middle', ...
            'fontname', 'courier', 'fontsize', 12);
end
drawnow;

function local_stdplot(x,y,dy,col, lw);
x = x(:); y = y(:); dy = dy(:);
dy(isnan(dy)) = 0;
vcol = (col+4*[1 1 1])/5;
patch([x; flipud(x)], [y+dy; flipud(y-dy)], vcol ,'edgecol', 'none', 'visible', 'off');
xplot(x,y,'color', col, 'linewidth' ,lw, 'visible', 'off');

function local_ssqplot(hs, AC);
Ncond = size(AC,1);
AC = AC(:);
axes(hs);LegStr = {};
for idif = 0:Ncond-1,
    ACd = AC([AC.jcond]-[AC.icond]==idif);
    hl = xplot([ACd.warp1], [ACd.totDiff], ['-' ploma(idif+1)], lico(idif+1));
    if idif==0, set(hl, 'linewidth', 2); end
    dwarp = ACd(1).warp2-ACd(1).warp1;
    LegStr = [LegStr, ['\Deltawarp=' num2str(dwarp) ' Oct.']];
end
hl = legend(LegStr, 'location', 'eastOutside', 'fontsize', 6);
apos = get(hs, 'position');
set(hl, 'units', 'normalized', 'position', [apos(1)+0.16 apos(2) 0.1 apos(4)]);
xlim([min([AC.warp1]) max([AC.warp1])]);
xlabel('Warp1');
ylabel('SSQ diff');
ylim([0 max(ylim)]);

function local_peakplot(hs, AC);
Ncond = size(AC,1);
AC = AC(:);
axes(hs);LegStr = {};
ACd = AC([AC.jcond]==[AC.icond]);
plot([ACd.warp1], [ACd.peak1Diff], '-ks', 'linewidth', 2, 'markersize',8);
xplot([ACd.warp1], [ACd.peak2Diff], '-ko', 'linewidth', 1);
xlabel('Warp (Oct.)');
ylabel('Peak(N+T)-Peak(N)');

function local_SACspec(D, BinWidth);
SPT = spiketimes(D, 1, 'no-unwarp');
totSAC = 0;
for icond=1:D.Stim.Presentation.Ncond,
    twarpfac = 2.^(-D.Stim.Warp(icond)); % duration factor
    AW = [0 D.Stim.ToneOnset*twarpfac];
    spt = AnWin(SPT(icond,:), AW);
    [SaC, Tau] = SPTCORR(spt,'nodiag', 40, BinWidth);
    totSAC = totSAC + SaC;
end
totSAC = totSAC(:) - mean(totSAC);
N = numel(totSAC);
dt = diff(Tau(1:2)); % sample period in ms
SACspec = A2dB(abs(fft(hann(N).*totSAC, 7*N)));
SACspec = SACspec(1:round(end/2));
df = 1e3/(8*N*dt); % freq spacing in Hz; factor eight from zero-padded FFT
dplot(df, SACspec);
ylim(10*ceil(max(SACspec/10))+[-40 0]);
[Fpeak, MgPeak, isort] = localmax(df, SACspec);
isort = isort(1:min(4,end)); % four largest peaks
[Fpeak,  MgPeak] = deal(Fpeak(isort),  MgPeak(isort));
PLM = 'ovsd';
for ii=1:numel(Fpeak),
    fr = Fpeak(ii);
    mg = MgPeak(ii);
    hl(ii) = xplot(fr,mg, PLM(ii), 'color', 'w', 'markerfacecolor' , 0.85*(ii-1)/4*[1 1 1]);
    LegStr{ii} = [num2str(round(fr)) ' Hz'];
    %text(fr, mg+2, [num2str(round(fr)) ' Hz']);
end
legend(hl, LegStr, 'location', 'northeast' ,'fontsize', 6)
