function S = JL_TKisi(Exp, Pen, ifreq, iSPL);
% JL_TKisi - AP interval analysis of TK data
%    Usage:
%      JL_TKisi(Exp, Pen, ifreq, iSPL);
%
%    Example: 
%        JL_TKisi('RG09181', 'pen2u3', 3, 1)
%
% To expedite, convert v7 mat files t pv6 ones: 
% dd = rdir('D:\USR\Thomas\data\deristat_results'); 
% RE = regexp({dd.fullname},'.*mat$', 'match'); 
% ie = cellfun(@isempty,RE);
% RE(ie)=[]; RE = [RE{:}];
% clear pen*
% for ii=1:numel(RE),
%   clear pen*
%   qq = load(RE{ii});
%   FdN = fieldnames(qq);
%     for jj=1:numel(FdN),
%       cmd = [FdN{jj} ' =  qq.' FdN{jj} ';']
%       eval(cmd);
%     end
%   save(RE{ii}, FdN{:}, '-v6');
% end
%
%  EXAMPLES
%     JL_TKisi('RG09174', 'pen2u3', 3, 1)

Rdir = 'D:\USR\Thomas\Data\deristat_results';
dd = dir(fullfile(Rdir, Exp, Pen, [Pen '*.mat']));
FFN = fullfile(Rdir, Exp, Pen, dd.name);
qq = load(FFN);
Fld = strrep(dd.name, '.mat', '');
qq = qq.(Fld);

freq = qq.carfreq(ifreq); % car freq in Hz
T = 1e3/freq;  % stim period in ms
SPL = qq.yval{iSPL};  % level in dB SPL
repdur = qq.repdur{iSPL};
if numel(repdur)>1,
    repdur = repdur(ifreq);
end

SPT = qq.spiketimes{iSPL, ifreq};
ISI = diff(SPT);
ISI(ISI>40) = nan;
SPT = rem(SPT, repdur);
SPT(SPT>qq.stimdur{iSPL}+5) = nan; % analysis window: stim + 5 ms after 
SPT(SPT<10) = nan; % analysis window starts 10 ms after stim onset
cSPT = mod(SPT,T);

set(figure,'units', 'normalized', 'position', [0.0922 0.487 0.341 0.41])
subplot(2,1,1);
plot(cSPT(2:end), ISI, 'ko', 'markersize',5);
xplot([0 T], [T T], 'r');
xplot([0 T], 2*[T T], 'r:');
xlim([0 T]);
set(gca,'fontsize',12);
ylog125([min(ISI)/1.5 50]);
ylabel('ISI (ms)');
title(sprintf('%s/%s: %d Hz,  %d dB SPL', Exp, Pen, round(freq), SPL), 'fontsize', 14);

subplot(2,1,2);
Nbin = min(round(numel(cSPT)/3), 41);
hist(cSPT, linspace(0,T,Nbin));
xlim([0 T]);
set(gca,'fontsize',12);
xlabel('Time re stim cycle (ms)');
ylabel('# Spikes per bin');


set(figure,'units', 'normalized', 'position', [0.52 0.347 0.373 0.408]);
plot(cSPT(1:end-1),cSPT(2:end),'d', 'markersize',5);
xlim([0 T])
ylim([0 T])
%set(gca,'fontsize',12);
xlabel('SPT(n) re stim cycle (ms)');
ylabel('SPT(n+1) re stim cycle (ms)');
grid on;

Tpeak = ginput;
Tpeak = Tpeak(:,1);
phPeak = exp(2*pi*i*Tpeak/T);
phPeak = phPeak(:).';
Npeak = numel(phPeak);

phSPT = exp(2*pi*i*cSPT/T);
phSPT = denan(phSPT);
phSPT = phSPT(:);

[phPeak, phSPT] = SameSize(phPeak, phSPT);
Dev = abs(phPeak - phSPT);
[dum, peak_ID] = min(Dev,[],2);

for ipeak=1:Npeak,
    for jpeak=1:Npeak,
        T(ipeak,jpeak) = sum(peak_ID(1:end-1)==ipeak & peak_ID(2:end)==jpeak);
    end
    N(ipeak,1) = sum(peak_ID==ipeak);
end
Tn = T/sum(T(:));
Ntot = sum(N);
Tmodel = N*N.'/Ntot;
S = CollectInStruct(T,N,Tn,Ntot, Tmodel);
T
Tmodel

