function XcycleHist(RW, AP, EPSP, AnWin, Nbin);
%  XcycleHist - cycle histogram based on AP/EPSP output of deristats
%     XcycleHist(RW, AP, EPSP, AnWin, Nbin);

if nargin<4, AnWin = []; end
if nargin<5, Nbin = 31; end
% get SGSR dataset
DS = sgsrdataset(RW.ExpID, RW.RecID);
SP = XgetStimprops(RW);

if isempty(AnWin),
    AnWin = [0 SP.burstDur];
end

tAP = [AP.tpeak];
tEPSP = [EPSP.tpeak, AP.tpeakEPSP]; % both failing EPSP and EPSPs embedded in APs

% relate time to stimulus onset
tAPsinceOnset = denan(tAP - mostRecent(tAP, SP.onset));
tEPSPsinceOnset = denan(tEPSP - mostRecent(tEPSP, SP.onset));
% select events within analysis window
tAPsinceOnset = tAPsinceOnset(betwixt(tAPsinceOnset, AnWin));
tEPSPsinceOnset = tEPSPsinceOnset(betwixt(tEPSPsinceOnset, AnWin));

Tstim = 1e3/SP.Fcar; % stimulus period in ms
phAP = rem(tAPsinceOnset, Tstim)/Tstim;
phEPSP = rem(tEPSPsinceOnset, Tstim)/Tstim;

Edges = linspace(0,1+1e-10,Nbin+1);
Nap = histc(phAP, Edges);
Nep = histc(phEPSP, Edges);

set(figure,'units', 'normalized', 'position', [0.28 0.346 0.58 0.525])
subplot(2,1,1);
h=bar(Edges, Nep, 'histc'); set(h, 'FaceColor', 'r'); hold on
h=bar(Edges+1, Nep, 'histc'); set(h, 'FaceColor', 'r')
title([DS.title '--' num2str(SP.Xval) ' ' SP.Xunit]);
xlim([0 2])
ylabel('# Events');
[Rep, PhiEPSP] = vectorstrength(tEPSPsinceOnset,SP.Fcar/1e3);
text(0.1, max(ylim)-0.1*diff(ylim), ['R=' num2str(Rep,2)])

subplot(2,1,2);
h=bar(Edges, Nap, 'histc'); hold on
h=bar(Edges+1, Nap, 'histc');
xlim([0 2])
xlabel('Phase (cycle)');
ylabel('# APs');
[Rap, PhiAP] = vectorstrength(tAPsinceOnset,SP.Fcar/1e3);
Dphi = mod(PhiAP-PhiEPSP,1);
Delay = Dphi*Tstim;
text(0.1, max(ylim)-0.1*diff(ylim), ['R=' num2str(Rap,2) '  \Deltat = ' num2str(Delay,3) ' ms'])





