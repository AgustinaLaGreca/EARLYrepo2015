function Y = JLspontPeaks2(Uidx, DTsmooth, MinDist_ms, ThrPeak)
% JLspontPeaks2 - counting and analyzing peaks of spontaneous activity
%   Y=JLspontPeaks2(Uidx, DTsmooth, MinDist_ms)

Fhighpass = 100; % Hz highpass filtering
MinDT = 0.3; % ms min dist  btw EPSPs

[DTsmooth, MinDist_ms, ThrPeak] = arginDefaults('DTsmooth/MinDist_ms/ThrPeak', 0.08, [], []);

MinDist_ms = replaceMatch(MinDist_ms, [], 0.5); % default 0.5-ms minimum distance
ThrPeak = replaceMatch(ThrPeak, [], 0.4); % mV default threshold peak size
% AP criterium
%APcrit = JL_APcrit(S);
S = JLgetRec(Uidx);
APcrit = [S.APthrSlope S.APwindow_start S.APwindow_end S.CutoutThrSlope];
if numel(APcrit)<4, APcrit(4) = APcrit(1); end  % default trunc Thr is AP detection thr

%read data, truncate spikes, highpass filter & smoothen
[D, DS] = readTKABF(S);
dt = D.dt_ms;
rec = D.AD(1,1).samples;
rec = rec-mean(rec(:));
[rec, tAP] = APtruncate2(rec, APcrit(4), dt, APcrit(2:3));
[B,A] = butter(5, 2e-3*Fhighpass*dt ,'high');
rec = filter(B,A,rec);
rec = smoothen(rec, DTsmooth, D.dt_ms);
if isequal('BFS', DS.stimtype),
    rec = rec(1:round(500/dt)); % first 500 ms
end

% find peaks that are at least MinDT ms apart
[Tpeak, Vpeak] = localmax(dt, rec, MinDT);
iok = Tpeak>1.1;

[Tpeak, Vpeak] = deal(Tpeak(iok), Vpeak(iok));

[dum Tsn Sn_left] = cutFromWaveform(dt, rec, Tpeak, [-MinDist_ms 0]);
[dum Tsn Sn_rite] = cutFromWaveform(dt, rec, Tpeak, [0  MinDist_ms]);
minLeft = min(Sn_left,[],1);
minRite = min(Sn_rite,[],1);
dsize(Vpeak, minLeft, ThrPeak)
ibig = find(Vpeak(:)'-minLeft>ThrPeak & Vpeak(:)'-minRite>ThrPeak);

% reduce peak collection to the big ones; extract snippets 
[Tpeak, Vpeak, minLeft, minRite] = deal(Tpeak(ibig), Vpeak(ibig), minLeft(ibig), minRite(ibig));
[dum Tsn Sn] = cutFromWaveform(dt, rec, Tpeak, [-0.5  0.5]);
Sn = Sn - SameSize(max(Sn,[],1), Sn);
MeanSnip = mean(Sn,2);
Tmid = 0.5*dt*numel(MeanSnip);
[dum dum narrowMeanSnip] = cutFromWaveform(dt, MeanSnip, Tmid, MinDist_ms*[-1 1]);
%deconvrec = deconv(rec,narrowMeanSnip);

[AC, ilag] = xcorr(rec, rec, 2000);
T_ac = ilag*dt;
%========PLOT===========
TTT = [num2str(D.icond) ': ' num2str(DS.xval(D.icond)) ' ' DS.x.Unit '  ' DS.filename '' DS.info];

fh1 = gcf;
set(fh1,'units', 'normalized', 'position', [0.0188 0.489 0.954 0.41]);
dplot(dt, rec); TracePlotInterface(gcf); 
%xdplot([dt MinDist_ms], deconvrec, 'g')
xplot(Tpeak, Vpeak, 'g*')
xplot(Tpeak, minLeft, 'mo')
xplot(Tpeak, minRite, 'm+')
fenceplot(Tpeak,ylim,'r');
title(TTT);

set(figure,'units', 'normalized', 'position', [0.583 0.161 0.257 0.777]);
plot(Tsn, Sn);
xplot(Tsn, MeanSnip, 'linewidth', 4, 'color', 0.7*[1 1 1]);

[ExpID, RecID, icond] = deal(D.ExpID, D.RecID, D.icond);

Y = CollectInStruct(ExpID, RecID, icond, TTT, '-', DTsmooth, MinDist_ms, ThrPeak, '-', ...
    dt, rec, T_ac, AC);





