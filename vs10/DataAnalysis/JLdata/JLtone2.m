function Vtail = JLtone2(ExpID, RecID, icond, t0, plotArg);
% JLtone2 - quick & dirty analysis of JL's FS data
%    JLtone(ExpID, RecID, icond, t0, plotArg);

if nargin<2, RecID=[]; end
if nargin<3, icond=[]; end
if nargin<4, t0=0; end
if nargin<5, plotArg=''; end

Tign = 20; % first Tign ms ignored

if numel(icond)>1, % combine multiple plots in one graph using recursion
    [icond, t0] = SameSize(icond, t0);
    for ii=1:numel(icond),
        Vtail{ii} = JLtone2(ExpID, RecID, icond(ii), t0(ii), ploco(ii));
    end
    return;
end

%---------single condition from here--------------

[D,DS]=readTKABF(ExpID, RecID, icond); 
Freq = round(DS.xval(icond));
T = 1e3/Freq; % ms stim period.
Dur = DS.BurstDur(1); % ms tone dur
SPL = DS.SPL;

Ncycle = floor((Dur-Tign)/T); % # stim cycles over which to average
AnDur = T*Ncycle; % analysis duration
AnWin = [Dur-AnDur Dur];

V=[D.AD(2:end,1).samples]; % recorded potential
dt = D.dt_ms;
tt = timeaxis(V(:,1), dt); % time axis in ms

Vtail = mean(V,2); % mean across reps
Vtail = Vtail(betwixt(tt,AnWin)); % V restricted to analysis window
% resample to enable average over stim cycles
NsamTail = numel(Vtail);
NsamTailNew = Ncycle*round(NsamTail/Ncycle); % # samples rounded to multiple of Ncycle
Vtail = interp1(linspace(0,1,NsamTail).', Vtail, linspace(0,1,NsamTailNew).');
Vtail = LoopMean(Vtail, -Ncycle);
VtailRep = repmat(Vtail, [Ncycle 1]);
dt_tail = dt*NsamTail/NsamTailNew;

rho = shufcorr(V);
Vmean = mean(V,2);
FreqStr = [num2str(Freq) ' Hz'];
if isequal('', plotArg),
    figure;
    xdplot([dt t0], V);
    xdplot([dt t0],Vmean, 'k', 'linewidth', 3);
    setGUIdata(gcf, 'LegStr', {FreqStr});
else,
    Lstr = cellify(getGUIdata(gcf, 'LegStr', {}));
    Lstr = [Lstr FreqStr];
    setGUIdata(gcf, 'LegStr', Lstr);
    %xdplot([dt t0],Vmean, plotArg);
    xdplot([dt_tail AnWin(1)+t0], VtailRep, plotArg);
end

Lstr = getGUIdata(gcf, 'LegStr', {});
legend(Lstr);
set(gcf,'units', 'normalized', 'position', [0.0344 0.516 0.902 0.384]);
TracePlotInterface(gcf); ylim auto; ylim(ylim);
title([DS.title '  '  num2str(Freq) ' Hz   ' trimspace(num2str(SPL))  ' dB  ; rho = ' num2str(rho,3)]);


