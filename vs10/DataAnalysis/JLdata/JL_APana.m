function [S, T] = JL_APana(Uidx, doPlot);
% JL_APana - waveform analysis of APs in JL data (EPSP-AP latency, etc)
%    [S, T] = JL_APana(Uidx, doPlot) 
%    S contains full detail of snippets, time derivatives, etc. T contains
%    NNTPs, which  are added to databes in file
%    'D:\processed_data\JL\NNTP\JL_APana.dbase'

[doPlot] = arginDefaults('doPlot',1);

tau_sm = 0.031; % ms smoothing window prior to diff
max_lat = 0.501; % ms max EPSP-AP latency
EPSPrisedur = 0.2; % ms max rising phase of EPSP precding EPSP peak or inflex
%minInfLat = 0.1;

CFN = fullfile('D:\processed_data\JL\AP_ANA_cache', [mfilename '_' num2str(Uidx)]);
[S_T, CFN, CP] = getcache(CFN, {Uidx tau_sm max_lat EPSPrisedur});
if ~isempty(S_T),
    [S, T] = deal(S_T.S, S_T.T);
else, % compute it
    [S, T] = local_compute(Uidx, tau_sm, max_lat, EPSPrisedur);
    S_T = CollectInStruct(S,T);
    % store in cache & database
    putcache(CFN, 1, CP, S_T);
    DBFN = 'D:\processed_data\JL\NNTP\JL_APana.dbase';
    if ~exist(DBFN, 'file'),
        init_dbase(DBFN, 'UniqueRecordingIndex');
    end
    Add2dbase(DBFN, T);
end

if doPlot,
    localPlot(S,T);
end

%================================================

function [S, T] = local_compute(Uidx, tau_sm, max_lat, EPSPrisedur);
W = JLwaveforms(Uidx);
SN = W.Snip(:,W.APinStim);
dSN = diff(smoothen(W.Snip, tau_sm, -W.dt))/W.dt;
d2SN = diff(smoothen(dSN, tau_sm, -W.dt))/W.dt;
tSN = W.Tsnip; 
dt = W.dt;
tdSN = W.Tsnip(2:end);

Nsp = size(SN,2);
qWin = tdSN<-0.03 & tdSN>-max_lat; % EPSPwindow as log index in tdSN
Inflex = false(1,Nsp);
[Lat, iEPSP, Vbaseline, Vepsp, dVepsp, Vpeak, Vcrit, maxEslope] = deal([]); % defaults for the no-AP case
%[Vbaseline, Vepsp, dVepsp, Vpeak, maxEslope] = 

for isp=1:Nsp,
    sn = SN(2:end, isp);
    dsn = diff(smoothen(SN(:, isp), tau_sm, -W.dt))/W.dt;
    dsn_sh = dsn([2:end end]);
    Vbaseline(1,isp) = mean(sn(tSN<tSN(1)+0.5));
    Vpeak(1,isp) = max(sn(betwixt(tSN, 0.05*[-1 1]))); % AP peak V
    Vcrit(1,isp) = (Vbaseline(isp)+Vpeak(isp))/2;
    d2sn = [0; d2SN(:,isp)];
    ipeak = find(dsn>0 & dsn_sh<=0 & sn>Vcrit(isp) & qWin, 1, 'last' );
    if isempty(ipeak), % find minimum of 2nd der
        ipeak = find(dsn(2:end-1)>=0 & d2SN(1:end-1,isp)<0 & d2SN(2:end,isp)>0 & qWin(1:end-2), 1, 'last' );
        Inflex(isp) = true;
        %plot(tdSN, dsn); xplot(tdSN(ipeak), dsn(ipeak), '*'); ; pause;
    end
    % if isp==17, keyboard; end
    if isempty(ipeak), 
        iEPSP(1,isp) = nan;
        Lat(1,isp) = nan;
        Vepsp(1,isp) = nan;
        dVepsp(1,isp) = nan;
        maxEslope(1,isp) = nan;
    else,
        iEPSP(1,isp) = ipeak;
        Lat(1,isp) = -tdSN(ipeak);
        Vepsp(1,isp) = sn(ipeak);
        Vbaseline(isp) = median(prctile(sn(1:ipeak-1), 10));
        dVepsp(1,isp) = dsn(ipeak);
        maxEslope(1,isp) = max(dsn(betwixt(tdSN, tdSN(ipeak)+[-EPSPrisedur 0])));
    end
%     if isp==64, %maxEslope(1,isp)<0,
%         figure(gcf); clf;
%         subplot(2,1,1);
%         plot(tSN(2:end), sn);
%         fenceplot(-Lat(1,isp), ylim, 'r');
%         title(num2str(isp));
%         subplot(2,1,2);
%         plot(tSN(2:end), dsn);
%         fenceplot(-Lat(1,isp), ylim, 'r');
%         keyboard;
%     end
end

[corr_EspspSlope, poly_EspspSlope,  slope_EspspSlope, Nlat] = local_CorrSlope(maxEslope, Lat);
[corr_rawVesp, poly_rawVesp, slope_rawVesp] = local_CorrSlope(Vepsp, Lat);
[corr_compVesp, poly_compVesp, slope_compVesp] = local_CorrSlope(Vepsp-Vbaseline, Lat);

UniqueRecordingIndex = Uidx;
S = CollectInStruct(UniqueRecordingIndex, '-EPSP_AP', dt, tSN, SN, dSN, d2SN, ...
    '-params', max_lat, tau_sm, EPSPrisedur, ...
    '-Latresults', Lat, iEPSP, Inflex, Vbaseline, Vepsp, dVepsp, Vpeak, maxEslope, Vcrit);


% mean stats
MeanLat = nanmean(Lat);
StdLat = nanstd(Lat);
Nspikes = numel(Lat);
T = CollectInStruct(UniqueRecordingIndex, '-EPSP_AP_stats', Nlat, MeanLat, StdLat, ...
    '-LatCorr', corr_EspspSlope, poly_EspspSlope,  slope_EspspSlope, ...
    corr_rawVesp, poly_rawVesp, slope_rawVesp, ...
    corr_compVesp, poly_compVesp, slope_compVesp);


%===========================
function [Corr, Poly, Slope, Nlat] = local_CorrSlope(X,Y);
iok = ~isnan(Y);
[X,Y] = deal(X(iok), Y(iok));
    if numel(X)<3,
    [Corr, Poly, Slope] = deal(nan);
    Nlat = 0;
else,
    Corr = corr(X(:), Y(:));
    Poly = polyfit(X,Y,1);
    Slope = Poly(1);
    Nlat = numel(X);
end

function localPlot(S,T);
set(figure,'units', 'normalized', 'position', [0.25 0.23 0.468 0.669]);
subplot(3,1,1);
plot((S.maxEslope)', (S.Lat)', '.');
h=xplot(xlim, polyval(T.poly_EspspSlope, xlim), 'r');
legend(h,sprintf('R=%1.2f', T.corr_EspspSlope));
ylim([0 S.max_lat]);
xlabel('max EPSP slope (V/s)');
ylabel('EPSP-AP Latency (ms)');
title(JLtitle(S.UniqueRecordingIndex));
%
subplot(3,1,2);
plot((S.Vepsp)', (S.Lat)', '.');
h=xplot(xlim, polyval(T.poly_rawVesp, xlim), 'r');
legend(h,sprintf('R=%1.2f', T.corr_rawVesp));
ylim([0 S.max_lat]);
xlabel('raw V(EPSP) (mV)')
ylabel('EPSP-AP Latency (ms)');
%
subplot(3,1,3);
plot((S.Vepsp-S.Vbaseline)', (S.Lat)', '.');
h=xplot(xlim, polyval(T.poly_compVesp, xlim), 'r');
legend(h,sprintf('R=%1.2f', T.corr_compVesp));
ylim([0 S.max_lat]);
xlabel('V(EPSP)-baseline (mV)')
ylabel('EPSP-AP Latency (ms)');















