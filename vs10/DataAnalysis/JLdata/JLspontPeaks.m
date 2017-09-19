function Y = JLspontPeaks(S, DTsmooth, PeakWidth, AcceptCrit, doPlot)
% JLspontPeaks - counting peaks of spontaneous activity
%   Y=JLspontPeaks(S, DTsmooth, PeakWidth, AcceptCrit, doPlot)

Fhighpass = 100; % Hz highpass filtering

if nargin<2, DTsmooth = []; end 
if nargin<3, PeakWidth = []; end 
if nargin<4, AcceptCrit = []; end
if nargin<5, doPlot = (nargout<1); end

if isempty(DTsmooth), DTsmooth = 0.2; end % ms default smoothing time
if isempty(PeakWidth), PeakWidth = 0.15; end % ms default smoothing time
if isempty(AcceptCrit), AcceptCrit = 0.3; end

[VspRaw, dt, tAP] = JLspont(S, 100, 0);
NspontAP = numel(tAP);
Vsp = smoothen(VspRaw, DTsmooth, dt);
[Tpeak, Vpeak] = localmax(dt, Vsp);
%
DVsp = smoothen(diff(Vsp,1)/dt, DTsmooth/2, dt);
NsamPeak = max(1,round(PeakWidth/dt));
DVspMax = runmax(DVsp, NsamPeak);
DVspMin = -runmax(-DVsp, NsamPeak);

D2Vsp = smoothen(diff(Vsp,2)/dt^2, DTsmooth/2, dt);
[Tpeak2, D2Vpeak] = localmax(dt, -D2Vsp);
D2Vpeak = -D2Vpeak;
ND2peak = numel(Tpeak2);


isamPeak = 1+round(Tpeak/dt);
DVmax = DVspMax(isamPeak);
DVmin = DVspMin(isamPeak);
DVtop = max(DVmax,-DVmin);
% reject smooth peaks
TopThr = AcceptCrit*median(DVtop);
iok = find(DVtop>TopThr);
[Tpeak, Vpeak, DVmax, DVmin, DVtop] = deal(Tpeak(iok), Vpeak(iok), DVmax(iok), DVmin(iok), DVtop(iok));
Npeak = numel(Tpeak);
SpontEventRate = 1e3*Npeak/(dt*numel(Vsp)); % #events/s

[ExpID, RecID, icond] = deal(S.ExpID, S.RecID, S.icond)
Y = CollectInStruct(ExpID, RecID, icond, '-', ...
    DTsmooth, PeakWidth, AcceptCrit, '-', ...
    tAP, NspontAP, '-', ...
    Tpeak, Vpeak, Npeak, SpontEventRate, '-', ...
    DVmax, DVmin, DVtop, '-', ...
    Tpeak2, D2Vpeak, ND2peak);


if doPlot,
    fh1 = figure;
    set(fh1,'units', 'normalized', 'position', [0.00625 0.534 0.987 0.376])
    JLfigmenu(fh1,S);
    fh2 = figure;
    set(fh2,'units', 'normalized', 'position', [0.00625 0.0859 0.987 0.376])
    JLfigmenu(fh2,S);

    figure(fh1);
    dplot(dt,VspRaw);
    xdplot(dt, Vsp, 'g');
    fenceplot(Tpeak,ylim,'r')

    figure(fh2);
    dplot([dt dt/2],DVsp);
    xdplot([dt dt/2],DVspMax,'m');
    xdplot([dt dt/2],DVspMin,'m');
    fenceplot(Tpeak,ylim,'r')
    xplot([0 500],[0 0], 'k');
    xplot(xlim, [1 1]*TopThr, 'color', 0.7*[1 1 1] ,'linewidth', 2);
    xplot(xlim, [-1 -1]*TopThr, 'color', 0.7*[1 1 1] ,'linewidth', 2);

    TracePlotInterface([fh1 fh2]);
    xlim([0 10])

%      fh3 = figure;
%      JLfigmenu(fh3,S);
%      set(fh3,'units', 'normalized', 'position', [0.822 0.337 0.17 0.237]);
%      hist(DVtop,50);
%      %      set(findobj(fh3,'type', 'patch'), 'facecolor', 'r')
%      %      hold on
%      %      hist(DVmin,50);
%      fh4 = figure;
%      plot(DVmax, -DVmin, '.');
%      xlim([0 max([xlim ylim] )]);
%      ylim([0 max(xlim)]);
%      set(gcf,'units', 'normalized', 'position', [0.541 0.261 0.262 0.298]);
%      xplot(xlim, xlim, 'k');
% 

    figure(fh1)
end



