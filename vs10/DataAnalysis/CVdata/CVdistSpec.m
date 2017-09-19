function CVdistSpec(DS, ichan, SPLshift)
% CVdistSpec - power spectrum of linear response cmps & 3rd-order DPs
%    CVdistSpec(DS, ichan, SPLshift)

[S, SP] = apple(DS, ichan);

f2i = @(f)1+round(f/SP.df); % freq [kHz] -> isample

Fprim = DS.Stim.Fzwuis/1e3; % stim freqs in kHz
DP = DPfreqs(Fprim, 3,1); % of 3-order, sumweight=1 DPs
DP = DP([DP.sumweight]==1 & [DP.order]==3); % reject sumweight -1 and primaries

Freq = Xaxis(SP.Cspec, SP.df); % freq axis in kHz

Magn = A2dB(abs(SP.Cspec)) + SPLshift;
DPfreq = cat(1,DP.freq);
isamDP = f2i(DPfreq);
isamPrim = f2i(Fprim);

DPmag = Magn(isamDP);
[DPfreq, DPmag] = sortAccord(DPfreq, DPmag, DPfreq);


% construct patch for outlining DP power
Xpatch = DPfreq([1 1:end end]);
Ypatch = [0; DPmag; 0];
ipatch = convhull(Xpatch, Ypatch);
[Xpatch, Ypatch] = deal(Xpatch(ipatch), Ypatch(ipatch));
% 
%patch(Xpatch, Ypatch, 'r');
stem(Freq(isamPrim), Magn(isamPrim), 'filled');
hold on
plot(DPfreq, DPmag, 'r.');
xlim([0 30]);
Ymax = ceil(max(10*Magn(isamPrim)))/10;
ylim([0 Ymax])
xlabel('Frequency');
ylabel('Power (dB)');
legend({'Linear response' '3-order DPs' }, 'location', 'northwest')
title(sprintf('%s rec %d   %d dB per component', expname(DS), irec(DS), DS.Stim.SPL));










