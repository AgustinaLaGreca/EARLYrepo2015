% JL_NNTP
% notes on database
DB = JLdbase; % returns all candidate recs 
DB300Hz = (DB([DB.freq]==300)); % subselection example
DB195_3 = DB([DB.iexp]==195 & [DB.icell]==3); % othr example

% detailed info on recordings
JLdatastruct(DB(23).UniqueRecordingIndex)

% get raw data, 50-Hz removed
[S, rawrec] = JLgetRec(DB195_3(11).UniqueRecordingIndex);

% grand database
GBdir = fullfile(processed_datadir,  '\JL\NNTP');

% view raw waveforms with some extras
JLviewrec(DB195_3(23).UniqueRecordingIndex); 
for ii=22:1e4, JLviewrec(DB195_3(ii).UniqueRecordingIndex); waitfor(gcf); pause(0.1); end

% =====this loop fills the JLspikeThr cache. First rm any caches
for ii=1:numel(DB),
    uix = DB(ii).UniqueRecordingIndex;
    JLspikeThr(uix);
    pause(0.2);
    aa;
end
%======this loop stores info to be retrieved by JLcycleStats & JLwaveforms=====
load(fullfile(GBdir, 'globalStats'), 'T')
istart = 1;
for ii=istart:numel(DB),
    uix = DB(ii).UniqueRecordingIndex;
    JLwaveforms(uix ,'rmcache'); % make sure cached waveforms are removed
    JLwaveforms(uix); % this forces recomputation of waveforms and cyclestats
end
% =====loop for JLcyclecorrs. Needs JLwaveforms and JLcycleStats to work, so first run previous loop =====
for ii=1:numel(DB),
    uix = DB(ii).UniqueRecordingIndex;
    CC(ii) = JLbeatVar2(uix);
    save(fullfile(GBdir, 'cycleCorrs'), 'CC', '-V6');
    disp([num2str(ii) '  ' JLtitle(uix)])
end
%======this loop stores info to be retrieved by JLanovaStats & JLanovaDetails & JLphaselock=====
PL = [];
for ii=1:numel(DB),
    uix = DB(ii).UniqueRecordingIndex;
    [A, Astats(ii)] = JLanova2(uix);
    pl = JLspikes(uix);
    PL = [PL pl];
    disp([num2str(ii) '  ' JLtitle(uix)])
    save(fullfile(GBdir, 'Phaselock'), 'PL', '-V6');
    aa; drawnow;
end

%=====================


%========combine different metrix=========
%=====AP_ANA database
%==XXXXX fix restriction to CanBeUsed
NNTP = structMerge('UniqueRecordingIndex', JLdbase, JL_rec_select, JL_rec_select, JLcycleStats, JLcycleCorrs, JLanovaStats, JLspikeStats, ...
    retrieve_dbase('D:\processed_data\JL\NNTP\PhaseLock.dbase'));
NNTP = NNTP([NNTP.CanBeUsed]); % restrict based on quality of rec and AP-detectability
T_LAT = retrieve_dbase('D:\processed_data\JL\NNTP\JL_APana.dbase');
NNTP = structMerge('UniqueRecordingIndex', NNTP, T_LAT);

%=======stats of AP latencies (all binaural stim)=========
NNTP_B = NNTP([NNTP.chan]=='B');
subplot(3,1,1); hist([NNTP_B.APlagI], linspace(0,5,50)); xlim([0 4])
subplot(3,1,2); hist([NNTP_B.APlagC], linspace(0,5,50)); xlim([0 4])
subplot(3,1,3); hist([NNTP_B.APlagC]-[NNTP_B.APlagI], linspace(-5,5,100)); xlim([-4 4])
nanmedian([NNTP_B.APlagI])
% ans =
%     0.2000
nanmedian([NNTP_B.APlagC])
% ans =
%     0.1776
[nanmedian([NNTP_B.APlagC]-[NNTP_B.APlagI]) nanstd([NNTP_B.APlagC]-[NNTP_B.APlagI]) ]
% ans =
%    -0.0273    1.4079
% No apparent asymmetry, suggesting we're seeing the same mix of I & C
% inputs as does the axon hillock.

% ========TRIPLETS========
SV = retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase');

% =======check: JLviewrec(223041608) . AP thr seems to be inconsistent w AP
% waveforms.

%============residual variance vs spont var======
set(figure,'units', 'normalized', 'position', [0.28 0.316 0.551 0.583])
subplot(2,2,1);
CS = JLcycleStats;
plot([CS.VarSpont1], [CS.BeatResVar], 'k.');
xlim([0 1]); ylim(xlim); xplot(xlim,xlim, 'k');
xlabel('Spont var (mV^2)')
ylabel('Bin res var (mV^2)')
subplot(2,2,2);
hist([CS.BeatResVar] - [CS.VarSpont1], linspace(-0.5,0.5,200)); xlim(0.3*[-1 1])
title('Residual driven var minus spont var')
xlabel('Var difference (mV^2)')
subplot(2,2,3);
plot([CS.IpsiResVar], [CS.BeatResVar], '.');
xplot([CS.ContraResVar], [CS.BeatResVar], '.r');
xlim([0 1]); ylim(xlim); xplot(xlim,xlim, 'k');
xlabel('Mon res var (mV^2)')
ylabel('Bin res var (mV^2)')
legend({'Ipsi', 'Contra'}, 'location', 'northwest')
subplot(2,2,4);
plot([CS.IpsiResVar]+[CS.ContraResVar], [CS.BeatResVar], '.');
xlabel('I+C res var (mV^2)')
ylabel('Bin res var (mV^2)')
xlim([0 1]); ylim(xlim); xplot(xlim,xlim, 'k');


%=================================
% baseline RMS: histogram of population
for icell=1:45,
    nntp = NNTP([NNTP.icell_run]==icell); % restrict to this cell
    MeanSpont1(icell) = mean([nntp.MeanSpont1]);
    MaxSpont1(icell) = max([nntp.VarSpont1]);
    VarSpont1(icell) = mean([nntp.VarSpont1]);
end
hist(sqrt(MaxSpont1),100);

% baseline RMS: time courses
Legstr = {};
for icell=1:45,
    nntp = NNTP([NNTP.icell_run]==icell); % restrict to this cell
    xplot([nntp.MinutesSinceRecStart], sqrt([nntp.VarSpont1]), [ploma(icell) '-'], lico(icell));
    Legstr{icell} = ['cell ' num2str(icell)];
end
legend(Legstr);



