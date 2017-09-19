function B = JLbinint(Yi, Yc, Yb, doPlot);
% JLbinint - analyze bin. interaction from I/C/B data
%   JLbinint(Si, Sc, Sb), where Sx are either JLbeat outputs (slow) or
%   JLbeatVar outputs (fast).
%
%   See also JLbeatvar.

doPlot = arginDefaults('doPlot', nargout<1);

if nargin==1, % Bflat output
    B = JLbinint(JLdatastruct(Yi.Uidx_I) ,JLdatastruct(Yi.Uidx_C) ,JLdatastruct(Yi.Uidx_B), doPlot);
    return;
end

if numel(Yi)>1, % array: handle recursively
    B = [];
    for ii=1:numel(Yi),
        b = JLbinint(Yi(ii), Yc(ii), Yb(ii), doPlot);
        B = [B, b];
    end
    return;
end

%===================single Yi, etc, from here==========

% check cache
Cdir = fullfile(processed_datadir, 'JL', 'JLbinint');
[B, CFN, CP] = getcache(fullfile(Cdir, 'BinIntCache'), local_recindex([Yi Yc Yb]));
if isempty(B), % compute
    if isnumeric(Yi), Yi = JLdatastruct(Yi); end
    if isnumeric(Yc), Yc = JLdatastruct(Yc); end
    if isnumeric(Yb), Yb = JLdatastruct(Yb); end

    if ~isfield(Yi, 'Freq1'), Yi = JLdatastruct(Yi.UniqueRecordingIndex); end
    if ~isfield(Yc, 'Freq1'), Yc = JLdatastruct(Yc.UniqueRecordingIndex); end
    if ~isfield(Yb, 'Freq1'), Yb = JLdatastruct(Yb.UniqueRecordingIndex); end

    if numel(unique([Yi.Freq1 Yc.Freq1 Yb.Freq1]))>1,
        error('Non identical frequencies.');
    end

    if numel(unique([Yi.SPL1 Yc.SPL2 Yb.SPL1 Yb.SPL2]))>1,
        error('Non identical SPLs.');
    end
    if ~isfield(Yi, 'bin_MeanWave') || isempty(Yb.ipsi_MeanWave), % beatplot format; convert to beatVar format
        Yi = JLbeatVar(Yi);
        Yc = JLbeatVar(Yc);
        Yb = JLbeatVar(Yb);
    end
    [Rmax_ipsi, tau_ipsi]    = local_corr(Yb.ipsi_dt,   Yb.ipsi_MeanWave, Yi.ipsi_MeanWave);
    [Rmax_contra, tau_contra]= local_corr(Yb.contra_dt, Yb.contra_MeanWave, Yc.contra_MeanWave);

    ipsi_Sil_totVar = Yi.sil_totVar;
    contra_Sil_totVar = Yc.sil_totVar;
    bin_Sil_totVar = Yb.sil_totVar;
    B = CollectInStruct(Yi, Yc, Yb, '-', Rmax_ipsi, tau_ipsi, Rmax_contra, tau_contra, ...
        '-', ipsi_Sil_totVar, contra_Sil_totVar, bin_Sil_totVar);
    putcache(CFN ,1e4, CP, B);
end

if ~doPlot, return; end
disp(sprintf('%% rho=[%1.2f %1.2f]; tau=[%i %i] us', B.Rmax_ipsi ,B.Rmax_contra, B.tau_ipsi, B.tau_contra));

%======================================================================
ThickGray = struct('color', 0.7*[1 1 1], 'linewidth', 3);

fh1 = figure;
set(fh1,'units', 'normalized', 'position', [0.0391 0.471 0.658 0.456])
figname(fh1, [mfilename ' mean&var vs time'])


% var
h1 = subplot(2,3,1);
dplot(B.Yi.ipsi_dt, B.Yi.ipsi_VarWave, 'b');
xlim([0 B.Yb.Period1]);
xplot(xlim, [1 1]*B.Yi.SpontVar, ThickGray);
ylabel('Across-Cycle variance (mV^2)');
h2 = subplot(2,3,2);
dplot(B.Yb.ipsi_dt, B.Yb.ipsi_VarWave, 'b');
xdplot(B.Yb.contra_dt, B.Yb.contra_VarWave, 'r');
xplot(xlim, [1 1]*B.Yb.SpontVar, ThickGray);
h3 = subplot(2,3,3);
dplot(B.Yc.ipsi_dt, B.Yc.contra_VarWave, 'r');
xplot(xlim, [1 1]*B.Yc.SpontVar, ThickGray);
ylim(h1, [0 max([ylim(h1) ylim(h2) ylim(h3)])]);
linkaxes([h1 h2 h3], 'xy')
title(B.Yb.TTT);
% mean
h4 = subplot(2,3,4);
dplot(B.Yi.ipsi_dt, B.Yi.ipsi_MeanWave, 'b');
xlim([0 B.Yb.Period1]);
%xplot(xlim, [1 1]*Yi.SpontMean, ThickGray);
xlabel('Time (ms)');
ylabel('Cycle-averaged V_m (mV)');
h5 = subplot(2,3,5);
dplot(B.Yb.ipsi_dt, B.Yb.ipsi_MeanWave, 'b');
xdplot(B.Yb.contra_dt, B.Yb.contra_MeanWave, 'r');
%xplot(xlim, [1 1]*B.Yb.SpontMean, ThickGray);
xlabel('Time (ms)');
h6 = subplot(2,3,6);
dplot(B.Yc.ipsi_dt, B.Yc.contra_MeanWave, 'r');
%xplot(xlim, [1 1]*Yc.SpontMean, ThickGray);
YL = [min([ylim(h4) ylim(h5) ylim(h6)]) max([ylim(h4) ylim(h5) ylim(h6)])];
ylim(h4, [YL(1) mean(YL)+0.8*diff(YL)]);
linkaxes([h4 h5 h6], 'xy')
legend(h5,{sprintf('R=%1.2f; delay %i \\mus', B.Rmax_ipsi, B.tau_ipsi), ...
    sprintf('R=%1.2f, delay %i \\mus', B.Rmax_contra, B.tau_contra)}, ...
    'location', 'North', 'fontsize', 9);
xlabel('Time (ms)');

fh2 = figure;
set(fh2,'units', 'normalized', 'position', [0.686 0.342 0.294 0.353]);
figname(fh2, [mfilename ' monomean  vs binmean'])
plot(B.Yi.ipsi_MeanWave, B.Yb.ipsi_MeanWave, 'b');
xplot(B.Yc.contra_MeanWave, B.Yb.contra_MeanWave, 'r');
LL = [min([xlim ylim]) max([xlim ylim])];
xlim(LL); ylim(LL); axis square; grid on
set(gca,'ytick', get(gca, 'xtick'));
xlabel('Monaural V_m (mV)');
ylabel('Binaural V_m (mV)');

% B.Yi = local_strip(B.Yi);
% B.Yc = local_strip(B.Yc);
% B.Yb = local_strip(B.Yb);


%==============================
function [rhomax, tau]=local_corr(dt, BinW, MonW);
Nmax = numel(BinW);
[rhomax, itau] = maxcorr(repmat(BinW, 10,1), repmat(MonW, 10,1), Nmax);
tau = round(1e3*dt*itau); % in mus

function S = local_strip(S);
for ii=1:numel(S),
    s = S(ii);
    s.bin_MeanWave = [];
    s.bin_VarWave = [];
    s.ipsi_MeanWave = [];
    s.ipsi_VarWave = [];
    s.contra_MeanWave = [];
    s.contra_VarWave = [];
    s.sil_MeanWave = [];
    s.sil_VarWave = [];
    S(ii) = s;
end


function Yr = local_recindex(Y);
for ii=1:numel(Y),
    if isstruct(Y(ii)),
        Yr(ii) = Y(ii).UniqueRecordingIndex;
    else,
        Yr(ii) = Y(ii);
    end
end















