function T = JLevents(Rec, DT, EPslopeThr, AnaWin, Ncat);
% JLevents - elementary event analysis of JLbeat data
%   JLEvents(Rec, DT, EPslopeThr, AnaWin, Ncat);
%         Rec: recording ID (see JLdatastruct)
%          DT: smoothing (ms) applied to dV/dt to estimate slopes [0.1 ms]
%  EPslopeThr: threshold EPSPS slope (V/s)
 
EventDur = 0.5; % ms typical event dur


[DT, EPslopeThr, AnaWin, Ncat] = arginDefaults('DT/EPslopeThr/AnaWin/Ncat', 0.1, 3, [0 inf], 12); 

% retrieve JLbeat output
S=JLdatastruct(Rec);

[D, DS, L]=readTKABF(S);
dt = D.dt_ms; % sample period in ms
NsamTot = D.NsamPerSweep/D.Nchan;
cutWin = [S.APwindow_start S.APwindow_end];

rec = D.AD(1).samples; % entire recording, including pre- and post-stim parts
drec = smoothen(diff(rec)/dt, DT, dt); % time derivative in V/s
rec = (rec(1:end-1)+rec(2:end))/2; % shift rec 1/2 sample to align with derivative
Time = Xaxis(rec, dt);
% restrict analyis to anawin =====
qana = betwixt(Time, AnaWin);
Time = Time(qana);
rec = rec(qana);
drec = drec(qana);
% [rec, SPTraw, Tsnip, Snip, APslope] = ...
%     APtruncate2(rec, [S.APthrSlope, S.CutoutThrSlope], dt, cutWin, 1); % last 1: replace APs by nans

% ===find events based on their upward slope===
NsamEvent = round(EventDur/dt);
drec_rm = runmax(drec, NsamEvent);
drec_rm(drec_rm<EPslopeThr) = min(drec); % discard low slopes
figure; hist(drec_rm, 200); xlim([0 inf]);
ievent = find(drec==drec_rm);
Tpeak_slope = Time(ievent);
% throw out events following the previous events too soon
itoosoon = 1+find(diff(Tpeak_slope)<EventDur);
itoosoon(find(diff(itoosoon)==1)+1) = []; % correct domino effect "short+short+short = short", etc
ievent(itoosoon) = [];
Tpeak_slope = Time(ievent);
Peak_slope = drec(ievent);
[dum, Tev, Events] = cutFromWaveform([dt, AnaWin(1)], rec, Tpeak_slope, [-1 2]*EventDur);
Events = Events - SameSize(mean(Events(betwixt(Tev,[-0.5 0.5]*EventDur),:),1), Events);
t2isam = @(t)1+round((t-AnaWin(1))/dt);
% ===categorize events based on the rising slope
I = categorize(Peak_slope, Ncat);

set(figure,'units', 'normalized', 'position', [0.0844 0.489 0.903 0.41]);
subplot(2,1,1)
plot(Time, rec);
xplot(Tpeak_slope, rec(ievent), 'r*')
subplot(2,1,2)
plot(Time, drec);
xplot(Time, drec_rm, 'g');
xplot(Tpeak_slope, Peak_slope, 'r*')
TracePlotInterface(gcf);

figure;
[]
CatPlot(Tev, Events, I);

figure;
isamev = [];
[rt0,rt1]= minmax(Tev)
for ii=1:numel(ievent),
    t_ev = Tpeak_slope(ii);
    i0 = t2isam(t_ev+rt0);
    i1 = t2isam(t_ev+rt1);
    if i0>0 && i1<numel(rec),
        isamev = [isamev, (i0:i1).'];
    end
end
isamev(isamev<1)=[];
isamev(isamev>numel(rec))=[];
xplot(rec(isamev),drec(isamev));






