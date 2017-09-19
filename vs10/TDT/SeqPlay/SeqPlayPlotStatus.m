function SeqPlayPlotStatus;
% SeqPlayPlotStatus - plot sequence play history (debug function)
%    SeqPlayPlotStatus plots the statistic of a seqplay action.
%    Usage:
%      ... (prepare sequenced play)
%      SeqPlayPlotStatus; % performs a seq

seqplaygo;
for ii=1:5e4,
    S(ii) = seqplaystatus;
    if ~S(ii).Active, break; end
end;
Nsam = length(S);
jitter = 0.1; 

set(gcf,'units', 'normalized', 'position', [0.475 0.12917 0.48563 0.7775])
subplot(3,1,1);
plot(cat(1,S.isam_abs), cat(1,S.isam_rel))

WaveCount = cat(1,S.WaveCount)+jitter*rand(Nsam,2);
WaveIndex = cat(1,S.WaveIndex)+jitter*rand(Nsam,2);
IplayList = cat(1,S.IplayList)+jitter*rand(Nsam,2);
Irep = cat(1,S.Irep)+jitter*rand(Nsam,2);

subplot(3, 1, 2);
plot(cat(1,S.isam_abs), WaveCount(:,1), 'b');
xplot(cat(1,S.isam_abs), WaveIndex(:,1), 'r');
xplot(cat(1,S.isam_abs), IplayList(:,1), 'g');
xplot(cat(1,S.isam_abs), Irep(:,1), 'm');
legend({'WaveCount'  'WaveIndex'  'IplayList' 'Irep'});


subplot(3, 1, 3);
plot(cat(1,S.isam_abs), WaveCount(:,2), 'b');
xplot(cat(1,S.isam_abs), WaveIndex(:,2), 'r');
xplot(cat(1,S.isam_abs), IplayList(:,2), 'g');
xplot(cat(1,S.isam_abs), Irep(:,2), 'm');





