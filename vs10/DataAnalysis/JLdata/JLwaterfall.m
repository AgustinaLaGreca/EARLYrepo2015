function JLwaterfall(S);
% JLwaterfall - waterfall plot of cyclic-averaged recordings
%   JLwaterfall(S), where S is struct array output of JLbeat, aligns the 
%   averaged recorded waveforms of different frequencies to their maximum
%   values. 

fh = figure;
set(fh,'units', 'normalized', 'position', [0.324 0.403 0.647 0.512])
[ah1, PhaseShift1, Tau1] = localplot(S, 1);
[ah2, PhaseShift2, Tau2] = localplot(S, 2);
YL1 = get(ah2,'ylim');
YL2 = get(ah2,'ylim');
YL = [min([YL1 YL2])  max([YL1 YL2])];
ylim(ah1,YL); ylim(ah2,YL);
xlabel('Time (ms)');
TracePlotInterface(gcf);
xlim([0 20]);

fh2 = figure; 
set(fh2,'units', 'normalized', 'position', [0.731 0.121 0.234 0.228]);
plot([S.Freq1], PhaseShift1, 'x-');
xplot([S.Freq2], PhaseShift2, 'ro-');
legend({['ipsi ' num2str(Tau1) ' ms'] ['contra ' num2str(Tau2) ' ms'] });
xlabel('Frequency (Hz)');
ylabel('Phase (cycle)');
figure(fh);

%fh3 = figure;


% =================
function [ah, PhaseShift, Tau] = localplot(S, ichan);
switch ichan,
    case 1, fn = 'Y1'; pn = 'Period1'; frn = 'Freq1'; Ttl ='ipsi';
    case 2, fn = 'Y2'; pn = 'Period2'; frn = 'Freq2'; Ttl ='contra';
end
ah = subplot(2,1,ichan);
for icond=1:numel(S),
    s = S(icond);
    Nsam(icond) = numel(s.Y1);
    tt = timeaxis(s.(fn), s.dt);
    Tmax = 2*max([s.Period1 s.Period2]); % max time plotted
    Y = s.(fn);
    Period = s.(pn);
    isam1 = tt<Period;
    [dum, imax] = max(Y(isam1));
    Tshift = s.dt*(imax-1) - Tmax/2 + 2*Period;
    plm = pmask(mod(tt,2*Period)<Period);
    xplot(tt-Tshift, Y, lico(icond));
    lh(icond) = xplot(tt-Tshift, Y+plm, lico(icond), 'linewidth', 1.5);
    LegStr{icond} = [num2str(s.(frn)) ' Hz'];
    isamshift(icond) = round(Tshift/s.dt);
    Tbeat = 1e3/s.Fbeat; % beat period
    PhaseShift(icond) = (s.dt*(imax-1)-0.5*Tbeat)/s.(pn);
end
Nsam = min(Nsam);
MeanY = [];
for icond=1:numel(S),
    y = circshift(S(icond).(fn), -isamshift(icond));
    MeanY = [MeanY, y(1:Nsam)];
end
MeanY = mean(MeanY,2);
lh(end+1) = xdplot(s.dt, MeanY, 'k--', 'linewidth',3);
LegStr{end+1} = 'mean';
legend(lh, LegStr);
title(Ttl, 'fontsize', 12);
xlim([0 Tmax]+0.25*Tmax);
freq = [S.(frn)]/1e3; % stim freq in kHz
[dum, Tau, PhaseShift] = LinPhaseFit(freq, PhaseShift, [2 10]);
Tau = deciRound(Tau,3);

