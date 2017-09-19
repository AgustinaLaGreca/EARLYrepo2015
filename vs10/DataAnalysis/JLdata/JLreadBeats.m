function [S, slist] = JLreadBeats(GerbilID, iCell, doCompute);
% JLreadBeats - read (or process) processed JLbeat data.
%   JLreadBeats(GerbilID, iCell, doCompute)
%   Processed data are load to Caller workspace.
%
%   JLreadBeats('rebuild') refreshes all caches using current settings.
%
%   [L, SL] = JLreadBeats('list') returns a master list of all JLbeat-processed 
%   MSO cells. SL is the same list as struct.
%


%====tools for determining AP thr============
% [CH RW AP EPSP] = deristats_tk(readTKaBF('RG10223', '3-1-BFS', 5), 0.05, -2);
% plot([AP.Vpeak]-[AP.Vdip],[AP.maxDownSlope], '.')
% or 
% mDslope = runmax(-RW.getdsam(), round(0.7/RW.dt)); hist(mDslope,100); ylim([0 3000])
%============================================

if nargin==1 && isequal('list', GerbilID), % rebuild cache
    S.RG09178 = [1 3]; 
    S.RG10191 = [1 2 3]; 
    S.RG10195 = 3; % 1,5 have only sinks, no positive going APs
    S.RG10196 = [2 4]; 
    S.RG10197 = [5 8]; 
    S.RG10198 = [1 2]; 
    S.RG10201 = [1 2 3]; 
    
    S.RG10204 = 10:13; 
    S.RG10209 = 2:6; % 7 has bookkeeping problems. Might be rescueable.
    S.RG10214 = 1:6;
    S.RG10216B = 1:2;
    S.RG10219 = 1:6;
    S.RG10223 = 1:7;
    slist = local_slist(S);
    return;
elseif nargin==1 && isequal('rebuild', GerbilID), % rebuild cache
    L = JLreadBeats('list');
    qqGerbilID = fieldnames(L);
    for ii=1:numel(qqGerbilID),
        ii
        gbid = qqGerbilID{ii}
        JLreadBeats(gbid, L.(gbid),1);
    end
    return;
end

if nargin<3,
    doCompute = 0;
end

if numel(iCell)>1, % call cells one by one using recursion
    for icell=iCell(:).';
        disp(['+++++ cell ' num2str(icell) ' +++++' ]);
        JLreadBeats(GerbilID, icell, doCompute);
    end
    return;
end

if isnumeric(GerbilID),
    GerbilID = gerbilName(GerbilID);
end

[SFN, ID] = local_filename(GerbilID, iCell); % filename for processed data
if doCompute,
    feval(['local_' ID], SFN);
end
xpr = ['load(''' SFN ''', ''Jb*'')'];
evalin('caller', 'clear Jb*');
evalin('caller', ['GerbilID = ''' GerbilID ''';']);
evalin('caller', ['iCell = ' num2str(iCell) ';']);
evalin('caller', xpr);


%===================================
function [SFN, ID] = local_filename(ExpID, iCell);
% filename for processed data
ID = [upper(ExpID) '_' num2str(iCell)];
SFN = fullfile(processed_datadir, 'JL', 'JLbeat', ID); % filename for processed data

%=============================
function local_addComments(cJL, cMH);
qq=evalin('caller', 'who(''Jb*'');');
for ii=1:numel(qq),
    cmd = ['[' qq{ii} '.varName] = deal(''' qq{ii} ''');'];
    evalin('caller', cmd);
    cmd = ['[' qq{ii} '.JLcomment] = deal(''' cJL ''');'];
    evalin('caller', cmd);
    cmd = ['[' qq{ii} '.MHcomment] = deal(''' cMH ''');'];
    evalin('caller', cmd);
end

%=======================================================================
%=======================================================================
%=======================================================================
%=======================================================================
%=======================================================================
%=======================================================================
%=======================================================================
%=======================================================================

function local_RG09178_1(SFN);
% RG10178 cell 1 ("EPSPs plus Aps");  Weakly binaural. Spikes only at 100 Hz. 
%      70 dB 2xbin
%-6    1-6-BFS  100*:150*:2500 Hz 70 dB     4 Hz beat                   17 x 1 x 3000/3500 ms   B  
%-38   1-38-BFS 100*:150*:2500 Hz 70 dB     4 Hz beat                   17 x 1 x 3000/3500 ms   B  
% many FS at different levels (!)
% 1 x CFS; % NTD, rho=1,-1,0; % BN
APcrit = [-15 -1 1.4 -12]; Nav = 1; 
for ifreq=1:17, 
    JbB_70dB(ifreq) = JLbeat('RG09178', '1-6-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
    JbB2_70dB(ifreq) = JLbeat('RG09178', '1-38-BFS', ifreq, Nav, APcrit, 0, 'B2_70dB');
end
local_addComments('EPSPs plus Aps', 'Weakly binaural. Spikes only at 100 Hz.');
save(SFN, 'Jb*'); clear Jb*;

function local_RG09178_3(SFN);
% RG10178 cell 3 ("big spikes, not much EPSPs?"); decent binaural cell
%      70 dB bin
%      70 dB, 100 Hz: beat freq varied
%      CFS, FS, NTD, everything @ 70 dB
% -54   3-7-BFS  100*:150*:2500 Hz 70 dB     4 Hz beat                   17 x 1 x 3000/3500 ms   B  
% -55   3-8-BFS  100|104 Hz        70 dB     4 Hz beat                   1 x 1 x 3000/3500 ms    B  
% -56   3-9-BFS  100|104 Hz        70 dB     4 Hz beat                   1 x 1 x 3000/3500 ms    B  
% -57   3-10-BFS 100|104 Hz        70 dB     4 Hz beat                   1 x 1 x 30000/35000 ms  B  
% -61   3-14-BFS 100|104 Hz        70 dB     4 Hz beat                   1 x 1 x 30000/35000 ms  B  
% -62   3-15-BFS 100|104 Hz        70 dB     4 Hz beat                   1 x 1 x 30000/35000 ms  B  
% -63   3-16-BFS 100 Hz            70 dB     0 Hz beat                   1 x 1 x 30000/35000 ms  B  
% -64   3-17-BFS 100|110 Hz        70 dB     10 Hz beat                  1 x 1 x 30000/35000 ms  B  
% -65   3-18-BFS 100|102 Hz        70 dB     2 Hz beat                   1 x 1 x 30000/35000 ms  B  
% -66   3-19-BFS 100 Hz            70 dB     0 Hz beat                   1 x 1 x 30000/35000 ms  B  
% -67   3-20-BFS 100|104 Hz        70 dB     4 Hz beat                   1 x 1 x 30000/35000 ms  B  
APcrit = [-20 -1 1.4 -15]; Nav = 1; % good crit
for ifreq=1:17, 
    JbB_70dB(ifreq) = JLbeat('RG09178', '3-7-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
end
local_addComments('big spikes, not much EPSPs?', 'decent binaural cell');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10191_1(SFN);
% RG10191 cell 1 ("big spikes, no clear EPSPs");  huge, atypical spikes. Only weakly binaural.
%      60 dB 2xbin/ipsi/contra
% XXX-1   1-1-BFS 100*:50*:600 Hz  60 dB     4 Hz beat 11 x 1 x 3000/3500 ms  B
% -2   1-2-BFS 100*:50*:600 Hz  60 dB     4 Hz beat 11 x 1 x 5000/5500 ms  B
% -3   1-3-BFS 100*:50*:600 Hz  60|-60 dB 4 Hz beat 11 x 1 x 5000/5500 ms  B
% -4   1-4-BFS 100*:50*:600 Hz  60|-60 dB 4 Hz beat 11 x 1 x 5000/6000 ms  B
% -5   1-5-BFS 100*:50*:600 Hz  -60|60 dB 4 Hz beat 11 x 1 x 5000/6000 ms  B
% -6   1-6-BFS 100*:50*:600 Hz  60 dB     4 Hz beat 11 x 1 x 5000/6000 ms  B
APcrit = [-12.5 -1 1.4]; Nav = 1; % good AP thr
for ifreq=1:10, 
    JbB_60dB(ifreq) = JLbeat('RG10191', '1-2-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
end
for ifreq=1:11,
    JbI_60dB(ifreq) = JLbeat('RG10191', '1-4-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10191', '1-5-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB2_60dB(ifreq) = JLbeat('RG10191', '1-6-BFS', ifreq, Nav, APcrit, 0, 'B2_60dB');
end
local_addComments('big spikes, no clear EPSPs', 'huge, atypical spikes. Only weakly binaural.');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10191_2(SFN);
% RG10191 cell 2 ("EPSP like data?")
%      60 dB 2xbin/ipsi/contra
%      CFS
% -9   2-3-BFS 100*:50*:1000 Hz 60 dB     4 Hz beat 19 x 1 x 9000/11000 ms B
% -10  2-4-BFS 100*:50*:1000 Hz -60|60 dB 4 Hz beat 19 x 1 x 9000/11000 ms B
% -11  2-5-BFS 100*:50*:1000 Hz 60|-60 dB 4 Hz beat 19 x 1 x 9000/11000 ms B
% XXX -12  2-6-BFS 100*:50*:1000 Hz 60 dB     4 Hz beat 19 x 1 x 9000/11000 ms B
APcrit = [-8 -1 1.4 -6]; Nav = 1; % ~good AP thr
for ifreq=1:19, 
    JbB_60dB(ifreq) = JLbeat('RG10191', '2-3-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10191', '2-4-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10191', '2-5-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    %JbB2_60dB(ifreq) = JLbeat('RG10191', '2-6-BFS', ifreq, Nav, APcrit, 0, '...');
end
local_addComments('EPSP like data?', '');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10191_3(SFN);
% RG10191 cell 3 ("1.5 mV spikes, not sure about EPSPs")
%      60 dB bin/ipsi/contra
%      CFS, FS, 
% -14  3-2-BFS 100*:50*:1000 Hz 60 dB     4 Hz beat 19 x 1 x 9000/11000 ms B
% -15  3-3-BFS 100*:50*:700 Hz  60|-60 dB 4 Hz beat 13 x 1 x 9000/11000 ms B
% -16  3-4-BFS 100*:50*:700 Hz  -60|60 dB 4 Hz beat 13 x 1 x 9000/11000 ms B
APcrit = [-8 -1 1.4]; Nav = 1; % good AP thr
for ifreq=1:19, 
    JbB_60dB(ifreq) = JLbeat('RG10191', '3-2-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
end
for ifreq=1:13, 
    JbI_60dB(ifreq) = JLbeat('RG10191', '3-3-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10191', '3-4-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
end
local_addComments('1.5 mV spikes, not sure about EPSPs', 'good AP thr');
save(SFN, 'Jb*'); clear Jb*;



function local_RG10195_3(SFN);
% RG10195 cell 3 ("EPSPs plus Aps or maybe multi unit?")  
%      65 dB bin/ipsi/contra
% FS, FM, BN
% -43  3-5-BFS  0.1*:0.1*:2 kHz  65 dB     4 Hz beat  20 x 1 x 6000/8500 ms  B  
% -44  3-6-BFS  0.1*:0.1*:1 kHz  75 dB     4 Hz beat  10 x 1 x 6000/8500 ms  B  
% -45  3-7-BFS  0.1*:0.1*:1 kHz  -75|75 dB 4 Hz beat  10 x 1 x 6000/8500 ms  B  
% -46  3-8-BFS  0.1*:0.1*:1 kHz  75|-75 dB 4 Hz beat  10 x 1 x 6000/8500 ms  B  
% -47  3-9-BFS  0.1*:0.1*:1 kHz  65|-65 dB 4 Hz beat  10 x 1 x 6000/8500 ms  B  
% -48  3-10-BFS 0.1*:0.1*:1 kHz  -65|65 dB 4 Hz beat  10 x 1 x 6000/8500 ms  B  
APcrit = [-10 -1 1.4]; Nav = 1; % well defined thr
for ifreq=1:20,
    JbB_65dB(ifreq) = JLbeat('RG10195','3-5-BFS', ifreq, Nav, APcrit, 0, 'B_65dB');
end
for ifreq=1:10,
    JbB_75dB(ifreq) = JLbeat('RG10195','3-6-BFS', ifreq, Nav, APcrit, 0, 'B_75dB');
    JbC_75dB(ifreq) = JLbeat('RG10195','3-7-BFS', ifreq, Nav, APcrit, 0, 'C_75dB');
    JbI_75dB(ifreq) = JLbeat('RG10195','3-8-BFS', ifreq, Nav, APcrit, 0, 'I_75dB');
    JbI_65dB(ifreq) = JLbeat('RG10195','3-9-BFS', ifreq, Nav, APcrit, 0, 'I_65dB');
    JbC_65dB(ifreq) = JLbeat('RG10195','3-10-BFS', ifreq, Nav, APcrit, 0, 'C_65dB');
end
local_addComments('EPSPs plus Aps or maybe multi unit?', 'well defined thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10196_2(SFN);
% RG10196 cell 2 ("EPSPs and Aps")  
%      80 dB bin/ipsi/contra
% FS
% -5   2-4-BFS 50*:50*:2000 Hz 80 dB     4 Hz beat 40 x 1 x 6000/8500 ms  B
% -6   2-5-BFS 50*:50*:1000 Hz 80|-80 dB 4 Hz beat 20 x 1 x 6000/8500 ms  B
% -7   2-6-BFS 50*:50*:1000 Hz -80|80 dB 4 Hz beat 20 x 1 x 6000/8500 ms  B
APcrit = [-10 -1 1.4 -7]; Nav = 1; % decent thr
for ifreq=1:40,
    JbB_80dB(ifreq) = JLbeat('RG10196','2-4-BFS', ifreq, Nav, APcrit, 0, 'B_80dB');
end
for ifreq=1:20,
    JbI_80dB(ifreq) = JLbeat('RG10196','2-5-BFS', ifreq, Nav, APcrit, 0, 'I_80dB');
    JbC_80dB(ifreq) = JLbeat('RG10196','2-6-BFS', ifreq, Nav, APcrit, 0, 'C_80dB');
end
local_addComments('EPSPs and Aps', 'decent thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10196_4(SFN);
% RG10196 cell 4 ("EPSPs and Aps")  
%      60 dB bin
%      70 dB bin: SINKS
%      80 dB bin/ipsi/contra
% FS
% -14  4-5-BFS 50*:50*:1000 Hz 80 dB     4 Hz beat 20 x 1 x 6000/8500 ms  B
% -15  4-6-BFS 50*:50*:1000 Hz -80|80 dB 4 Hz beat 20 x 1 x 6000/8500 ms  B
% -16  4-7-BFS 50*:50*:1000 Hz 80|-80 dB 4 Hz beat 20 x 1 x 6000/8500 ms  B
% -17  4-8-BFS 50*:50*:1000 Hz 60 dB     4 Hz beat 20 x 1 x 6000/8500 ms  B
% -18  4-9-BFS 50*:50*:1000 Hz 70 dB     4 Hz beat 20 x 1 x 6000/8500 ms  B
APcrit = [-10 -1 1.4 -7]; Nav = 1; % mediocre thr definition
for ifreq=1:20,
    JbB_80dB(ifreq) = JLbeat('RG10196','4-5-BFS', ifreq, Nav, APcrit, 0, 'B_80dB');
    JbC_80dB(ifreq) = JLbeat('RG10196','4-6-BFS', ifreq, Nav, APcrit, 0, 'C_80dB');
    JbI_80dB(ifreq) = JLbeat('RG10196','4-7-BFS', ifreq, Nav, APcrit, 0, 'I_80dB');
    JbB_60dB(ifreq) = JLbeat('RG10196','4-8-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbB_70dB(ifreq) = JLbeat('RG10196','4-9-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
end
local_addComments('EPSPs and Aps', 'mediocre thr definition');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10197_5(SFN);
% RG10197 cell 5 ("big spikes")  
%      70 dB 2xbin/ipsi/contra
% -5   5-1-BFS 50*:50*:1000 Hz 70 dB     4 Hz beat 20 x 1 x 6000/8500 ms  B
% -6   5-2-BFS 50*:50*:1000 Hz -70|70 dB 4 Hz beat 20 x 1 x 6000/8500 ms  B
% -7   5-3-BFS 50*:50*:1000 Hz 70|-70 dB 4 Hz beat 20 x 1 x 6000/8500 ms  B
% -8   5-4-BFS 50*:50*:200 Hz  70 dB     4 Hz beat 4 x 1 x 6000/8500 ms   B
APcrit = [-20 -1 1.4 -13]; Nav = 1; % pretty reliable thr definition
for ifreq=1:20,
    JbB_70dB(ifreq) = JLbeat('RG10197','5-1-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
    JbC_70dB(ifreq) = JLbeat('RG10197','5-2-BFS', ifreq, Nav, APcrit, 0, 'C_70dB');
    JbI_70dB(ifreq) = JLbeat('RG10197','5-3-BFS', ifreq, Nav, APcrit, 0, 'I_70dB');
end
for ifreq=1:4,
    JbB2_70dB(ifreq) = JLbeat('RG10197','5-4-BFS', ifreq, Nav, APcrit, 0, 'B2_70dB');
end
local_addComments('big spikes', 'pretty reliable thr definition');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10197_8(SFN);
% RG10197 cell 8 ("big spikes")  
%      70 dB 2xbin/ipsi/contra
%      FS
% -19  8-1-BFS 50*:50*:1000 Hz 70 dB     4 Hz beat 20 x 1 x 6000/8500 ms  B
% -20  8-2-BFS 100*:50*:600 Hz 70|-70 dB 4 Hz beat 11 x 1 x 6000/8500 ms  B
% -21  8-3-BFS 100*:50*:600 Hz -70|70 dB 4 Hz beat 11 x 1 x 6000/8500 ms  B
% -22  8-4-BFS 100*:50*:600 Hz 70 dB     4 Hz beat 11 x 1 x 6000/8500 ms  B
APcrit = [-20 -1 1.4]; Nav = 1; % immense spikes; trivial to trigger them.
for ifreq=1:20,
    JbB_70dB(ifreq) = JLbeat('RG10197','8-1-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
end
for ifreq=1:11,
    JbI_70dB(ifreq) = JLbeat('RG10197','8-2-BFS', ifreq, Nav, APcrit, 0, 'I_70dB');
    JbC_70dB(ifreq) = JLbeat('RG10197','8-3-BFS', ifreq, Nav, APcrit, 0, 'C_70dB');
    JbB2_70dB(ifreq) = JLbeat('RG10197','8-4-BFS', ifreq, Nav, APcrit, 0, 'B2_70dB');
end
local_addComments('big spikes', 'immense spikes; trivial to trigger them.');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10198_1(SFN);
% RG10198 cell 1 ("spikes + EPSPs")  
%      70 dB bin/contra
% -1   1-2-BFS  0.1*:0.1*:2 kHz  70 dB     4 Hz beat               20 x 1 x 6000/8500 ms B
% -2   1-3-BFS  0.1*:0.1*:2 kHz  -70|70 dB 4 Hz beat               20 x 1 x 6000/8500 ms B
APcrit = [-10 -1 1.4]; Nav = 1; % good AP thr
for ifreq=1:11, % one missing abf screws up the bookkeeping
    JbB_70dB(ifreq) = JLbeat('RG10198','1-2-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
end
for ifreq=1:20,
    JbC_70dB(ifreq) = JLbeat('RG10198','1-3-BFS', ifreq, Nav, APcrit, 0, 'C_70dB');
end
local_addComments('spikes + EPSPs', 'good AP thr; bookkeeping error: JbB_70dB misses one ABF');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10198_2(SFN);
% RG10198 cell 2 ("spikes + EPSPs + heartbeat")  
%      40 dB 3xbin/ipsi/contra
%      70 dB 2xbin/ipsi/2xcontra
%      FS
% -3   2-1-BFS  0.1*:0.1*:2 kHz  70 dB     4 Hz beat               20 x 1 x 6000/8500 ms B
% -4   2-2-BFS  50*:50*:1000 Hz  70 dB     4 Hz beat               20 x 1 x 6000/8500 ms B
% -5   2-3-BFS  50*:50*:1000 Hz  70|-70 dB 4 Hz beat               20 x 1 x 6000/8500 ms B
% -6   2-4-BFS  50*:50*:1000 Hz  -70|70 dB 4 Hz beat               20 x 1 x 6000/8500 ms B
% -17  2-15-BFS 50*:50*:1000 Hz  40 dB     4 Hz beat               20 x 1 x 6000/8500 ms B
% -18  2-16-BFS 50*:50*:1000 Hz  -40|40 dB 4 Hz beat               20 x 1 x 6000/8500 ms B
% -19  2-17-BFS 50*:50*:1000 Hz  40|-40 dB 4 Hz beat               20 x 1 x 6000/8500 ms B
% -20  2-18-BFS 50*:50*:1000 Hz  40 dB     4 Hz beat               20 x 1 x 6000/8500 ms B
% -21  2-19-BFS 50*:50*:1000 Hz  40 dB     4 Hz beat               20 x 1 x 6000/8500 ms B
% -22  2-20-BFS 50*:50*:1000 Hz  -40|40 dB 4 Hz beat               20 x 1 x 6000/8500 ms B
APcrit = [-10 -1 1.4]; Nav = 1; % good AP thr
for ifreq=1:20, 
    JbB_70dB(ifreq) = JLbeat('RG10198','2-1-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
    JbB2_70dB(ifreq) = JLbeat('RG10198','2-2-BFS', ifreq, Nav, APcrit, 0, 'B2_70dB');
    JbI_70dB(ifreq) = JLbeat('RG10198','2-3-BFS', ifreq, Nav, APcrit, 0, 'I_70dB');
    JbC_70dB(ifreq) = JLbeat('RG10198','2-4-BFS', ifreq, Nav, APcrit, 0, 'C_70dB');
    JbB_40dB(ifreq) = JLbeat('RG10198','2-15-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbC_40dB(ifreq) = JLbeat('RG10198','2-16-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
    JbI_40dB(ifreq) = JLbeat('RG10198','2-17-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
    JbB2_40dB(ifreq) = JLbeat('RG10198','2-18-BFS', ifreq, Nav, APcrit, 0, 'B2_40dB');
    JbB3_40dB(ifreq) = JLbeat('RG10198','2-19-BFS', ifreq, Nav, APcrit, 0, 'B3_40dB');
    JbC2_40dB(ifreq) = JLbeat('RG10198','2-20-BFS', ifreq, Nav, APcrit, 0, 'C2_40dB');
end
local_addComments('spikes + EPSPs + heartbeat', 'good AP thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10201_1(SFN);
% RG10201 cell 1 ("spikes + EPSPs")  
%      60 dB bin/ipsi
% -1   1-1-BFS  250*:150*:2500 Hz 60 dB     4 Hz beat 16 x 1 x 6000/8500 ms B
% -2   1-2-BFS  250*:150*:2500 Hz 60|-60 dB 4 Hz beat 16 x 1 x 6000/8500 ms B
APcrit = [-14 -1 1.4]; Nav = 1; % good AP thr
for ifreq=1:16, 
    JbB_60dB(ifreq) = JLbeat('RG10201', '1-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10201', '1-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
end
local_addComments('spikes + EPSPs', 'good AP thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10201_2(SFN);
% RG10201 cell 2 ("spikes + EPSPs")  
%      40 dB 2xbin/ipsi/contra
%      60 dB bin [xt freq]
%      60 dB bin/ipsi/contra
% -3   2-1-BFS  250*:150*:2500 Hz 60 dB     4 Hz beat 16 x 1 x 6000/8500 ms B
% -4   2-2-BFS  50*:100*:1050 Hz  60 dB     4 Hz beat 11 x 1 x 6000/8500 ms B
% -5   2-3-BFS  50*:100*:1050 Hz  60|-60 dB 4 Hz beat 11 x 1 x 6000/8500 ms B
% -6   2-4-BFS  50*:100*:1050 Hz  -60|60 dB 4 Hz beat 11 x 1 x 6000/8500 ms B
% -7   2-5-BFS  50*:100*:1050 Hz  40 dB     4 Hz beat 11 x 1 x 6000/8500 ms B
% XXXXX -8   2-6-BFS  50*:100*:1050 Hz  -40|40 dB 4 Hz beat 2* x 1 x 6000/8500 ms B
% XXXXX -9   2-7-BFS  50*:100*:1050 Hz  -40|40 dB 4 Hz beat 2* x 1 x 6000/8500 ms B
% -10  2-8-BFS  50*:100*:1050 Hz  -40|40 dB 4 Hz beat 11 x 1 x 6000/8500 ms B
% -11  2-9-BFS  50*:100*:1050 Hz  40|-40 dB 4 Hz beat 11 x 1 x 6000/8500 ms B
% -12  2-10-BFS 50*:100*:1050 Hz  40 dB     4 Hz beat 11 x 1 x 6000/8500 ms B
APcrit = [-14.5 -1 1.4]; Nav = 1; % good AP thr
for ifreq=1:16, 
    JbB_60xdB(ifreq) = JLbeat('RG10201', '2-1-BFS', ifreq, Nav, APcrit, 0, 'B_60xdB');
end
for ifreq=1:11, 
    JbB_60dB(ifreq) = JLbeat('RG10201', '2-2-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10201', '2-3-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10201', '2-4-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB_40dB(ifreq) = JLbeat('RG10201', '2-5-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
%     JbC_40dB(ifreq) = JLbeat('RG10201', '2-6-BFS', ifreq, Nav, APcrit, 0, '...');
%     JbC2_40dB(ifreq) = JLbeat('RG10201', '2-7-BFS', ifreq, Nav, APcrit, 0, '...');
    JbC_40dB(ifreq) = JLbeat('RG10201', '2-8-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
    JbI_40dB(ifreq) = JLbeat('RG10201', '2-9-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
    JbB2_40dB(ifreq) = JLbeat('RG10201', '2-10-BFS', ifreq, Nav, APcrit, 0, 'B2_40dB');
end
local_addComments('spikes + EPSPs', 'good AP thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10201_3(SFN);
% RG10201 cell 3 ("spikes + EPSPs")  
%      40 dB bin/contra
%      60 dB bin
% -13  3-1-BFS  50*:100*:1050 Hz  60 dB     4 Hz beat 11 x 1 x 6000/8500 ms B
% -14  3-2-BFS  1.05:0.1:2.05 kHz 60 dB     4 Hz beat 11 x 1 x 6000/8500 ms B
% -15  3-3-BFS  50*:100*:2050 Hz  40 dB     4 Hz beat 21 x 1 x 6000/8500 ms B
% -16  3-4-BFS  50*:100*:2050 Hz  -40|40 dB 4 Hz beat 21 x 1 x 6000/8500 ms B
% -17  3-5-BFS  50*:100*:2050 Hz  40|-40 dB 4 Hz beat 6* x 1 x 6000/8500 ms B
APcrit = [-12.5 -1 1.4]; Nav = 1; % good AP thr
for ifreq=1:11, 
    JbB_60dB(ifreq) = JLbeat('RG10201', '3-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
end
for ifreq=2:11, 
    JbB_60dB(ifreq+10) = JLbeat('RG10201', '3-2-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
end
for ifreq=1:21, 
    JbB_40dB(ifreq) = JLbeat('RG10201', '3-3-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbC_40dB(ifreq) = JLbeat('RG10201', '3-4-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
end
local_addComments('spikes + EPSPs', 'good AP thr');
save(SFN, 'Jb*'); clear Jb*;


function local_RG10204_10(SFN);
% RG10204 cell 10 ("spikes + EPSPs")  
%      50 dB bin/ipsi/contra
%      60 dB 3xbin/2xipsi/2xcontra
%      70 dB bin/ipsi/contra
%      71 dB 6xbin/7xipsi/6xcontra (70 dB, 5 freqs, 300-700 Hz)
%      FS, NTD, BN, FM, CFS
% -57   10-1-BFS   100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -58   10-2-BFS   100*:100*:1500 Hz 60|-60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -59   10-3-BFS   100*:100*:1500 Hz -60|60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -60   10-4-BFS   100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -61   10-5-NTD   100~1000 Hz       60 dB     -1500:100:1500 us rho=1  31 x 10 x 200/300 ms   B  
% -62   10-6-NTD   100~1000 Hz       60 dB     -1500:100:1500 us rho=-1 31 x 10 x 200/300 ms   B  
% -63   10-7-NTD   100~1000 Hz       60 dB     -1500:100:1500 us rho=0  31 x 10 x 200/300 ms   B  
% -64   10-8-FS    75.4::9000 Hz     60 dB     no mod                   24 x 10 x 70/100 ms    B  
% -65   10-9-FS    75.4::9000 Hz     60|-60 dB no mod                   24 x 10 x 70/100 ms    B  
% -66   10-10-FS   75.4::9000 Hz     -60|60 dB no mod                   24 x 10 x 70/100 ms    B  
% -67   10-11-FS   75.4::9000 Hz     -60|40 dB no mod                   24 x 10 x 70/100 ms    B  
% -68   10-12-FS   75.4::9000 Hz     40|-40 dB no mod                   24 x 10 x 70/100 ms    B  
% -69   10-13-FS   75.4::9000 Hz     40 dB     no mod                   24 x 10 x 70/100 ms    B  
% -70   10-14-FS   75.4::9000 Hz     30 dB     no mod                   24 x 10 x 70/100 ms    B  
% -71   10-15-FS   75.4::9000 Hz     -30|30 dB no mod                   24 x 10 x 70/100 ms    B  
% -72   10-16-FS   75.4::9000 Hz     30|-30 dB no mod                   24 x 10 x 70/100 ms    B  
% -73   10-17-FS   75.4::9000 Hz     -50|50 dB no mod                   24 x 10 x 70/100 ms    B  
% -74   10-18-FS   75.4::9000 Hz     50 dB     no mod                   24 x 10 x 70/100 ms    B  
% -75   10-19-FS   75.4::9000 Hz     50|-50 dB no mod                   24 x 10 x 70/100 ms    B  
% -76   10-20-FS   75.4::9000 Hz     70 dB     no mod                   24 x 10 x 70/100 ms    B  
% -77   10-21-FS   75.4::9000 Hz     70 dB     no mod                   24 x 10 x 70/100 ms    B  
% -78   10-22-FS   75.4::9000 Hz     70|-70 dB no mod                   24 x 10 x 70/100 ms    B  
% -79   10-23-FS   75.4::9000 Hz     -70|70 dB no mod                   24 x 10 x 70/100 ms    B  
% -80   10-24-FM   *****             80 dB     ***                      1 x 1 x 10500/11000 ms B  
% -81   10-25-FM   *****             80 dB     ***                      1 x 1 x 10500/11000 ms B  
% -82   10-26-FM   *****             80 dB     ***                      1 x 1 x 10500/11000 ms B  
% -83   10-27-FM   *****             80 dB     ***                      1 x 1 x 10500/11000 ms B  
% -84   10-28-FM   *****             80 dB     ***                      1 x 1 x 10500/11000 ms B  
% -85   10-29-FM   *****             70 dB     ***                      1 x 1 x 10500/11000 ms B  
% -86   10-30-FM   *****             70 dB     ***                      1 x 1 x 10500/11000 ms B  
% -87   10-31-FM   *****             70 dB     ***                      1 x 1 x 10500/11000 ms B  
% -88   10-32-FM   *****             70 dB     ***                      1 x 1 x 10500/11000 ms B  
% -89   10-33-FM   *****             70 dB     ***                      1 x 1 x 10500/11000 ms B  
% -90   10-34-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -91   10-35-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -92   10-36-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -93   10-37-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -94   10-38-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -95   10-39-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -96   10-40-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -97   10-41-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -98   10-42-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -99   10-43-FM   *****             60 dB     ***                      1 x 1 x 10500/11000 ms B  
% -100  10-44-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -101  10-45-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -102  10-46-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -103  10-47-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -104  10-48-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -105  10-49-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -106  10-50-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -107  10-51-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -108  10-52-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -109  10-53-FM   *****             60|-60 dB ***                      1 x 1 x 10500/11000 ms B  
% -110  10-54-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -111  10-55-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -112  10-56-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -113  10-57-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -114  10-58-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -115  10-59-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -116  10-60-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -117  10-61-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -118  10-62-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -119  10-63-FM   *****             -60|60 dB ***                      1 x 1 x 10500/11000 ms B  
% -120  10-64-CFS  20 Hz             30 dB     500 us +                 1 x 51 x 60/100 ms     B  
% -121  10-65-CFS  20 Hz             40 dB     500 us +                 1 x 51 x 60/100 ms     B  
% -122  10-66-CFS  20 Hz             50 dB     500 us +                 1 x 51 x 60/100 ms     B  
% -123  10-67-CFS  20 Hz             60 dB     500 us +                 1 x 51 x 60/100 ms     B  
% -124  10-68-CFS  20 Hz             70 dB     500 us +                 1 x 51 x 60/100 ms     B  
% 1     10-69-BN   950 Hz            45 dB     20/90/0.25               1 x 1 x 79/79 s        B/L
% 2     10-70-BN   950 Hz            45 dB     20/90/0.25               1 x 1 x 64/64 s        B/R
% 3     10-71-BN   950 Hz            45 dB     20/90/0.25               1 x 1 x 79/79 s        B/L
% 4     10-72-BN   950 Hz            45 dB     20/90/0.25               1 x 1 x 64/64 s        B/B
% 5     10-73-BN   950 Hz            45 dB     20/90/0.25               1 x 1 x 64/64 s        B/B
% 6     10-74-BN   950 Hz            45 dB     20/90/0.25               1 x 1 x 64/64 s        B/R
% -125  10-75-BFS  100*:100*:1500 Hz 50 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -126  10-76-BFS  100*:100*:1500 Hz -50|50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -127  10-77-BFS  100*:100*:1500 Hz 50|-50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -128  10-78-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -129  10-79-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -130  10-80-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -131  10-81-BFS  100*:100*:1500 Hz 70 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -132  10-82-BFS  100*:100*:1500 Hz -70|70 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -133  10-83-BFS  100*:100*:1500 Hz 70|-70 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -134  10-84-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -135  10-85-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -136  10-86-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -137  10-87-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -138  10-88-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -139  10-89-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -140  10-90-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -141  10-91-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -142  10-92-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -143  10-93-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -144  10-94-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -145  10-95-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -146  10-96-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -147  10-97-BFS  300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -148  10-98-BFS  300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -149  10-99-BFS  300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -150  10-100-BFS 300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -151  10-101-BFS 300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -152  10-102-BFS 300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -153  10-103-NTD 100~3000 Hz       70 dB     0 us rho=0               1 x 100 x 200/300 ms   B  
% -154  10-104-NTD 100~3000 Hz       60 dB     0 us rho=0               1 x 100 x 200/300 ms   B  
% -155  10-105-NTD 100~3000 Hz       60|-60 dB 0 us rho=0               1 x 100 x 200/300 ms   B  
% -156  10-106-NTD 100~3000 Hz       -60|60 dB 0 us rho=0               1 x 100 x 200/300 ms   B  
% -157  10-107-NTD 100~3000 Hz       60|-60 dB 0 us rho=0               1 x 100 x 200/300 ms   B  

% -57   10-1-BFS   100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -58   10-2-BFS   100*:100*:1500 Hz 60|-60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -59   10-3-BFS   100*:100*:1500 Hz -60|60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -60   10-4-BFS   100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B APcrit = [-16 -1 1.4]; Nav = 1; % reasonably well-defined thr
% -125  10-75-BFS  100*:100*:1500 Hz 50 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -126  10-76-BFS  100*:100*:1500 Hz -50|50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -127  10-77-BFS  100*:100*:1500 Hz 50|-50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -128  10-78-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -129  10-79-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -130  10-80-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -131  10-81-BFS  100*:100*:1500 Hz 70 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -132  10-82-BFS  100*:100*:1500 Hz -70|70 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -133  10-83-BFS  100*:100*:1500 Hz 70|-70 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
%
% -134  10-84-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -135  10-85-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -136  10-86-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -137  10-87-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -138  10-88-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -139  10-89-BFS  300*:100*:700 Hz  70 dB     4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -140  10-90-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -141  10-91-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -142  10-92-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -143  10-93-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -144  10-94-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -145  10-95-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -146  10-96-BFS  300*:100*:700 Hz  -70|70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -147  10-97-BFS  300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -148  10-98-BFS  300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -149  10-99-BFS  300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -150  10-100-BFS 300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -151  10-101-BFS 300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
% -152  10-102-BFS 300*:100*:700 Hz  70|-70 dB 4 Hz beat                5 x 1 x 6000/10000 ms  B  
APcrit = [-10 -1 1.4]; Nav = 1; % AP thr is amazingly stable over recordings
for ifreq=1:15,
    JbB0_60dB(ifreq) = JLbeat('RG10204','10-1-BFS', ifreq, Nav, APcrit, 0, 'B0_60dB');% see JL's remark
    JbI_60dB(ifreq) = JLbeat('RG10204','10-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10204','10-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB_60dB(ifreq)= JLbeat('RG10204','10-4-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbB_50dB(ifreq) = JLbeat('RG10204','10-75-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10204','10-76-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10204','10-77-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbB2_60dB(ifreq) = JLbeat('RG10204','10-78-BFS', ifreq, Nav, APcrit, 0, 'B2_60dB');
    JbI2_60dB(ifreq) = JLbeat('RG10204','10-79-BFS', ifreq, Nav, APcrit, 0, 'I2_60dB');
    JbC2_60dB(ifreq) = JLbeat('RG10204','10-80-BFS', ifreq, Nav, APcrit, 0, 'C2_60dB');
    JbB_70dB(ifreq) = JLbeat('RG10204','10-81-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
    JbC_70dB(ifreq) = JLbeat('RG10204','10-82-BFS', ifreq, Nav, APcrit, 0, 'C_70dB');
    JbI_70dB(ifreq) = JLbeat('RG10204','10-83-BFS', ifreq, Nav, APcrit, 0, 'I_70dB');
end
for ifreq=1:5,
    JbB_71dB(ifreq) = JLbeat('RG10204','10-84-BFS', ifreq, Nav, APcrit, 0, 'B_71dB');
    JbB2_71dB(ifreq) = JLbeat('RG10204','10-85-BFS', ifreq, Nav, APcrit, 0, 'B2_71dB');
    JbB3_71dB(ifreq) = JLbeat('RG10204','10-86-BFS', ifreq, Nav, APcrit, 0, 'B3_71dB');
    JbB4_71dB(ifreq) = JLbeat('RG10204','10-87-BFS', ifreq, Nav, APcrit, 0, 'B4_71dB');
    JbB5_71dB(ifreq) = JLbeat('RG10204','10-88-BFS', ifreq, Nav, APcrit, 0, 'B5_71dB');
    JbB6_71dB(ifreq) = JLbeat('RG10204','10-89-BFS', ifreq, Nav, APcrit, 0, 'B6_71dB');
    JbC_71dB(ifreq) = JLbeat('RG10204','10-90-BFS', ifreq, Nav, APcrit, 0, 'C_71dB');
    JbC2_71dB(ifreq) = JLbeat('RG10204','10-91-BFS', ifreq, Nav, APcrit, 0, 'C2_71dB');
    JbC3_71dB(ifreq) = JLbeat('RG10204','10-92-BFS', ifreq, Nav, APcrit, 0, 'C3_71dB');
    JbC4_71dB(ifreq) = JLbeat('RG10204','10-93-BFS', ifreq, Nav, APcrit, 0, 'C4_71dB');
    JbC5_71dB(ifreq) = JLbeat('RG10204','10-94-BFS', ifreq, Nav, APcrit, 0, 'C5_71dB');
    JbC6_71dB(ifreq) = JLbeat('RG10204','10-95-BFS', ifreq, Nav, APcrit, 0, 'C6_71dB');
    JbC7_71dB(ifreq) = JLbeat('RG10204','10-96-BFS', ifreq, Nav, APcrit, 0, 'C7_71dB');
    JbI_71dB(ifreq) = JLbeat('RG10204','10-97-BFS', ifreq, Nav, APcrit, 0, 'I_71dB');
    JbI2_71dB(ifreq) = JLbeat('RG10204','10-98-BFS', ifreq, Nav, APcrit, 0, 'I2_71dB');
    JbI3_71dB(ifreq) = JLbeat('RG10204','10-99-BFS', ifreq, Nav, APcrit, 0, 'I3_71dB');
    JbI4_71dB(ifreq) = JLbeat('RG10204','10-100-BFS', ifreq, Nav, APcrit, 0, 'I4_71dB');
    JbI5_71dB(ifreq) = JLbeat('RG10204','10-101-BFS', ifreq, Nav, APcrit, 0, 'I5_71dB');
    JbI6_71dB(ifreq) = JLbeat('RG10204','10-102-BFS', ifreq, Nav, APcrit, 0, 'I6_71dB');
end
local_addComments('spikes + EPSPs', 'AP thr is amazingly stable over recordings');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10204_11(SFN);
% RG10204 cell 11 ("EPSPs + very small spikes")  
%      40 dB bin/ipsi/contra
%      50 dB bin/ipsi/contra
%      60 dB bin/ipsi/contra
% -158  11-1-BFS   100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -159  11-2-BFS   100*:100*:1500 Hz 60|-60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -160  11-3-BFS   100*:100*:1500 Hz -60|60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -161  11-4-BFS   100*:100*:1500 Hz -60|50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -162  11-5-BFS   100*:100*:1500 Hz 50|-50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -163  11-6-BFS   100*:100*:1500 Hz 50 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -164  11-7-BFS   100*:100*:1500 Hz 40 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -165  11-8-BFS   100*:100*:1500 Hz 40|-40 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -166  11-9-BFS   100*:100*:1500 Hz -40|40 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
APcrit = [-8 -1 1.4 -5.5]; Nav = 1; % poorly defined thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10204','11-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10204','11-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10204','11-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbC_50dB(ifreq) = JLbeat('RG10204','11-4-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10204','11-5-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbB_50dB(ifreq) = JLbeat('RG10204','11-6-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbB_40dB(ifreq) = JLbeat('RG10204','11-7-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbI_40dB(ifreq) = JLbeat('RG10204','11-8-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
    JbC_40dB(ifreq) = JLbeat('RG10204','11-9-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
end
local_addComments('EPSPs + very small spikes', 'poorly defined thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10204_12(SFN);
% RG10204 cell 12 ("spikes + EPSPs")  
%      30 dB bin/contra
%      40 dB bin/ipsi/contra
%      50 dB bin/ipsi/contra
%      60 dB bin/ipsi/contra
% -167  12-1-BFS   100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -168  12-2-BFS   100*:100*:1500 Hz 60|-60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -169  12-3-BFS   100*:100*:1500 Hz -60|60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -170  12-4-BFS   100*:100*:1500 Hz -60|50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -171  12-5-BFS   100*:100*:1500 Hz 50 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -172  12-6-BFS   100*:100*:1500 Hz 50|-50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -173  12-7-BFS   100*:100*:1500 Hz 40|-50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -174  12-8-BFS   100*:100*:1500 Hz 40 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -175  12-9-BFS   100*:100*:1500 Hz -40|40 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -176  12-10-BFS  100*:100*:1500 Hz -40|30 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -177  12-11-BFS  100*:100*:1500 Hz 30 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
APcrit = [-5 -1 1.4 -3.75]; Nav = 1; % reasonably well defined thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10204','12-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10204','12-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10204','12-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbC_50dB(ifreq) = JLbeat('RG10204','12-4-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbB_50dB(ifreq) = JLbeat('RG10204','12-5-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10204','12-6-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbI_40dB(ifreq) = JLbeat('RG10204','12-7-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
    JbB_40dB(ifreq) = JLbeat('RG10204','12-8-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbC_40dB(ifreq) = JLbeat('RG10204','12-9-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
    JbC_30dB(ifreq) = JLbeat('RG10204','12-10-BFS', ifreq, Nav, APcrit, 0, 'C_30dB');
    JbB_30dB(ifreq) = JLbeat('RG10204','12-11-BFS', ifreq, Nav, APcrit, 0, 'B_30dB');
end
local_addComments('spikes + EPSPs', 'reasonably well defined thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10204_13(SFN);
% RG10204 cell 13 ("spikes + EPSPs", "strong beat")  
%      40 dB bin
%      50 dB bin/ipsi/contra
%      60 dB bin/ipsi/contra
% NTD
% -178  13-1-BFS   100*:100*:1500 Hz 60 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -179  13-2-BFS   100*:100*:1500 Hz 60|-60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -180  13-3-BFS   100*:100*:1500 Hz -60|60 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -186  13-9-BFS   100*:100*:1500 Hz 50 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
% -187  13-10-BFS  100*:100*:1500 Hz 50|-50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -188  13-11-BFS  100*:100*:1500 Hz -50|50 dB 4 Hz beat                15 x 1 x 6000/10000 ms B  
% -189  13-12-BFS  100*:100*:1500 Hz 40 dB     4 Hz beat                15 x 1 x 6000/10000 ms B  
APcrit = [-9 -1 1.4 -7]; Nav = 1; % reasonably well defined thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10204','13-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10204','13-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10204','13-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB_50dB(ifreq) = JLbeat('RG10204','13-9-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10204','13-10-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10204','13-11-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbB_40dB(ifreq) = JLbeat('RG10204','13-12-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
end
local_addComments('spikes + EPSPs; strong beat', 'reasonably well defined thr');
save(SFN, 'Jb*'); clear Jb*;



function local_RG10209_2(SFN);
% RG10209 cell 2 ("spikes + EPSPs")  
%      60 dB bin/ipsi/contra
% -3   2-2-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B
% -4   2-3-BFS 100*:100*:1500 Hz 60|-60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
% -5   2-4-BFS 100*:100*:1500 Hz -60|60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
APcrit = [-20 -1 1.4 -15]; Nav = 1; % reasonably well-defined thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10209','2-2-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10209', '2-3-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10209', '2-4-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
end
local_addComments('spikes + EPSPs', 'reasonably well defined thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10209_3(SFN);
% RG10209 cell 3 ("spikes + EPSPs" "isolation changed during recording")  
%      60 dB 3xbin/ipsi/contra
% XXX-6   3-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B artifacts 
% -7   3-2-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B isolation changes during recording, better use 3-3-
% -8   3-3-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B signal still changing
% -9   3-4-BFS 100*:100*:1500 Hz -60|60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B signal more stable, but smaller
% -10  3-5-BFS 100*:100*:1500 Hz 60|-60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
APcrit = [-20 -1 1.4 -15]; Nav = 1; % reasonably well-defined thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10209','3-2-BFS', ifreq, Nav, [-35 -1 1.4 -20], 0, 'B_60dB');  %huge spikes
    JbB2_60dB(ifreq) = JLbeat('RG10209','3-3-BFS', ifreq, Nav, [-18 -1 1.4], 0, 'B2_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10209','3-4-BFS', ifreq, Nav, [-27 -1 1.4 -20], 0, 'C_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10209','3-5-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
end
local_addComments('spikes + EPSPs; isolation changed during recording', 'reasonably well-defined thr; JbB_60dB has huge spikes');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10209_4(SFN);
% RG10209 cell 4 ("spikes + EPSPs", "beats")  
%      60 dB 2xbin/ipsi/contra
% -11  4-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B
% -12  4-2-BFS 100*:100*:1500 Hz -60|60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
% -13  4-3-BFS 100*:100*:1500 Hz 60|-60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
% -14  4-4-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B
% 2 x NTD
APcrit = [-20 -1 1.4 -15]; Nav = 1; % reasonably well-defined thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10209', '4-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB'); 
    JbC_60dB(ifreq) = JLbeat('RG10209', '4-2-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10209', '4-3-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbB2_60dB(ifreq) = JLbeat('RG10209','4-4-BFS', ifreq, Nav, APcrit, 0, 'B2_60dB');
end
local_addComments('spikes + EPSPs; beats', 'reasonably well-defined thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10209_5(SFN);
% RG10209 cell 5 ("spikes + EPSPs, lost cell after 500 Hz")  
%      60 dB bin
% -17  5-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              7* x 1 x 6000/10000 ms B
APcrit = [-20 -1 1.4 -17]; Nav = 1; % reasonably well-defined thr
for ifreq=1:7,
    JbB_60dB(ifreq) = JLbeat('RG10209', '5-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB'); 
end
local_addComments('spikes + EPSPs, lost cell after 500 Hz', 'reasonably well-defined thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10209_6(SFN);
% RG10209 cell 6 ("lost cell")  
%      60 dB bin
% -18  6-2-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              3* x 1 x 6000/10000 ms B
APcrit = [-7.5 -1 1.4 -5]; Nav = 1; % reasonably well-defined thr
for ifreq=1:3,
    JbB_60dB(ifreq) = JLbeat('RG10209', '6-2-BFS', ifreq, Nav, APcrit, 0, 'B_60dB'); 
end
local_addComments('lost cell', 'reasonably well-defined thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10209_7(SFN);
% RG10209 cell 7 ("spikes + EPSPs", "instable signal")  
%      60 dB bin
% -19  7-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B
% -20  7-2-BFS 100*:100*:1500 Hz 60|-60 dB 4 Hz beat              9* x 1 x 6000/10000 ms B
% -21  7-3-BFS 100*:100*:1500 Hz -60|60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
error('Bookkeeping problems');
APcrit = [-7.5 -1 1.4]; Nav = 1; % reasonably well-defined thr
for ifreq=1:3,
    JbB_60dB(ifreq) = JLbeat('RG10209', '7-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB'); 
    JbI_60dB(ifreq) = JLbeat('RG10209', '7-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB'); 
    JbC_60dB(ifreq) = JLbeat('RG10209', '7-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB'); 
end
local_addComments('spikes + EPSPs; instable signal', 'Bookkeeping problems');
save(SFN, 'Jb*'); clear Jb*;



function local_RG10214_1(SFN);
% RG10214 cell 1 ("spikes + EPSPs" "lost cell during recording")  
%      60 dB bin
% -1   1-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
APcrit = [-6 -1 1.4 -5]; Nav = 1; % reasonably well-defined thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10214', '1-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
end
local_addComments('spikes + EPSPs; lost cell during recording', 'reasonably well-defined thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10214_2(SFN);
% RG10214 cell 2 ("strange waveform", "beats")  
%      60 dB 2xbin/ipsi/contra
% -2   2-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -3   2-2-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -4   2-3-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -5   2-4-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
APcrit = [-8 -1 1.4 -6]; Nav = 1; % hard to define a thr w this inverted (?) waveform
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10214', '2-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10214', '2-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10214', '2-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB2_60dB(ifreq) = JLbeat('RG10214', '2-4-BFS', ifreq, Nav, APcrit, 0, 'B2_60dB');
end
local_addComments('strange waveform; beats', 'hard to define a thr w this inverted (?) waveform');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10214_3(SFN);
% RG10214 cell 3 ("EPSPS no clear Aps")  
%      60 dB bin/ipsi/contra
% -6   3-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -7   3-2-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -8   3-3-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
APcrit = [-7.5 -1 1.4 -6]; Nav = 1; % not a sharp thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10214', '3-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10214', '3-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10214', '3-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
end
local_addComments('EPSPS no clear Aps', 'not a sharp thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10214_4(SFN);
% RG10214 cell 4 ("spikes + EPSPs", "beats and inhibition", "heart beats?")  
%      60 dB 2xbin/ipsi/contra
% -10  4-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -11  4-2-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -12  4-3-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -13  4-4-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
APcrit = [-7.5 -1 1.4 -6]; Nav = 1; % not a sharp thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10214',  '4-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10214',  '4-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10214',  '4-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB2_60dB(ifreq) = JLbeat('RG10214', '4-4-BFS', ifreq, Nav, APcrit, 0, 'B2_60dB');
end
local_addComments('spikes + EPSPs; beats and inhibition; heart beats', 'not a sharp thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10214_5(SFN);
% RG10214 cell 5 ("losing cell")  
%      60 dB bin
% -14  5-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
APcrit = [-4.5 -1 1.4 -3.5]; Nav = 1; % strange waveforms: sinks
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10214',  '5-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
end
local_addComments('losing cell', 'strange waveforms: sinks');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10214_6(SFN);
% RG10214 cell 6 ("spikes + EPSPs")  
%      10 dB bin/ipsi/contra
%      20 dB bin/ipsi/contra
%      30 dB bin/ipsi/contra
%      40 dB bin/ipsi/contra
%      50 dB bin/ipsi/contra
%      60 dB 3xbin/ipsi/contra
%      70 dB bin/ipsi/contra
%      80 dB bin/ipsi/contra
%   NTD
% -15  6-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -16  6-2-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -17  6-3-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -18  6-4-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -19  6-5-BFS  100*:100*:1500 Hz 70 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -20  6-6-BFS  100*:100*:1500 Hz -70|70 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -21  6-7-BFS  100*:100*:1500 Hz 70|-70 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -22  6-8-BFS  100*:100*:1500 Hz 80 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -23  6-9-BFS  100*:100*:1500 Hz -80|80 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -24  6-10-BFS 100*:100*:1500 Hz 80|-80 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -25  6-11-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -26  6-12-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -27  6-13-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -28  6-14-BFS 100*:100*:1500 Hz 40 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -29  6-15-BFS 100*:100*:1500 Hz 40|-40 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -30  6-16-BFS 100*:100*:1500 Hz -40|40 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -31  6-17-BFS 100*:100*:1500 Hz 30 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -32  6-18-BFS 100*:100*:1500 Hz -30|30 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -33  6-19-BFS 100*:100*:1500 Hz 30|-30 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -34  6-20-BFS 100*:100*:1500 Hz 20 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -35  6-21-BFS 100*:100*:1500 Hz 20|-20 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -36  6-22-BFS 100*:100*:1500 Hz -20|20 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -37  6-23-BFS 100*:100*:1500 Hz 10 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -38  6-24-BFS 100*:100*:1500 Hz 10|-10 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -39  6-25-BFS 100*:100*:1500 Hz -10|10 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -40  6-26-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
APcrit = [-12 -1 1.4 -10]; Nav = 1; % thr okayish
for ifreq=1:15,
    ifreq
    JbB_60dB(ifreq) = JLbeat('RG10214',  '6-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10214',  '6-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10214',  '6-3-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB2_60dB(ifreq) = JLbeat('RG10214', '6-4-BFS', ifreq, Nav, APcrit, 0, 'B2_60dB');
    JbB_70dB(ifreq) = JLbeat('RG10214',  '6-5-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
    JbC_70dB(ifreq) = JLbeat('RG10214',  '6-6-BFS', ifreq, Nav, APcrit, 0, 'C_70dB');
    JbI_70dB(ifreq) = JLbeat('RG10214',  '6-7-BFS', ifreq, Nav, APcrit, 0, 'I_70dB');
    JbB_80dB(ifreq) = JLbeat('RG10214',  '6-8-BFS', ifreq, Nav, APcrit, 0, 'B_80dB');
    JbC_80dB(ifreq) = JLbeat('RG10214',  '6-9-BFS', ifreq, Nav, APcrit, 0, 'C_80dB');
    JbI_80dB(ifreq) = JLbeat('RG10214',  '6-10-BFS', ifreq, Nav, APcrit, 0, 'I_80dB');
    JbB_50dB(ifreq) = JLbeat('RG10214',  '6-11-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10214',  '6-12-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10214',  '6-13-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbB_40dB(ifreq) = JLbeat('RG10214',  '6-14-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbI_40dB(ifreq) = JLbeat('RG10214',  '6-15-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
    JbC_40dB(ifreq) = JLbeat('RG10214',  '6-16-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
    JbB_30dB(ifreq) = JLbeat('RG10214',  '6-17-BFS', ifreq, Nav, APcrit, 0, 'B_30dB');
    JbC_30dB(ifreq) = JLbeat('RG10214',  '6-18-BFS', ifreq, Nav, APcrit, 0, 'C_30dB');
    JbI_30dB(ifreq) = JLbeat('RG10214',  '6-19-BFS', ifreq, Nav, APcrit, 0, 'I_30dB');
    JbB_20dB(ifreq) = JLbeat('RG10214',  '6-20-BFS', ifreq, Nav, APcrit, 0, 'B_20dB');
    JbI_20dB(ifreq) = JLbeat('RG10214',  '6-21-BFS', ifreq, Nav, APcrit, 0, 'I_20dB');
    JbC_20dB(ifreq) = JLbeat('RG10214',  '6-22-BFS', ifreq, Nav, APcrit, 0, 'C_20dB');
    JbB_10dB(ifreq) = JLbeat('RG10214',  '6-23-BFS', ifreq, Nav, APcrit, 0, 'B_10dB');
    JbI_10dB(ifreq) = JLbeat('RG10214',  '6-24-BFS', ifreq, Nav, APcrit, 0, 'I_10dB');
    JbC_10dB(ifreq) = JLbeat('RG10214',  '6-25-BFS', ifreq, Nav, APcrit, 0, 'C_10dB');
    JbB3_60dB(ifreq) = JLbeat('RG10214', '6-26-BFS', ifreq, Nav, APcrit, 0, 'B3_60dB');
end
local_addComments('spikes + EPSPs', 'okayish');
save(SFN, 'Jb*'); clear Jb*;


function local_RG10216B_1(SFN);
% RG216B cell 1 ("mini spikes, large EPSPS")  
%      70 dB 2xbin/contra/ipsi % few spikes, but good phase locking
% -1   1-1-BFS  100*:100*:1500 Hz 70 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
% -2   1-2-BFS  100*:100*:1500 Hz 70|-70 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -3   1-3-BFS  100*:100*:1500 Hz -70|70 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -4   1-4-BFS  100*:100*:1500 Hz 70 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
APcrit = [-7 -1 1 -5.5]; Nav = 1;% reasonable thr def
for ifreq=1:15,
    JbB_70dB(ifreq) = JLbeat('RG10216b', '1-1-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
    JbB2_70dB(ifreq) = JLbeat('RG10216b', '1-4-BFS', ifreq, Nav, APcrit, 0, 'B2_70dB');
    JbI_70dB(ifreq) = JLbeat('RG10216b', '1-2-BFS', ifreq, Nav, APcrit, 0, 'I_70dB');
    JbC_70dB(ifreq) = JLbeat('RG10216b', '1-3-BFS', ifreq, Nav, APcrit, 0, 'C_70dB');
end
local_addComments('mini spikes, large EPSPS', 'reasonable thr def');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10216B_2(SFN);
% RG216B cell 2 ("mini spikes, large EPSPS")
%      40 dB bin/contra/ipsi 
%      50 dB bin/contra/ipsi 
%      60 dB 5xbin/contra/ipsi 
%      70 dB bin/contra/ipsi 
% -5   2-1-BFS  100*:100*:1500 Hz 70 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
% -6   2-2-BFS  100*:100*:1500 Hz -70|70 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -7   2-3-BFS  100*:100*:1500 Hz 70|-70 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -8   2-4-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
% -9   2-5-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -10  2-6-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -11  2-7-BFS  100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
% -12  2-8-BFS  100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -13  2-9-BFS  100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -14  2-10-BFS 100*:100*:1500 Hz 40 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
% -15  2-11-BFS 100*:100*:1500 Hz 40|-40 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -16  2-12-BFS 100*:100*:1500 Hz -40|40 dB 4 Hz beat 15 x 1 x 6000/10000 ms B
% -17  2-13-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
% -18  2-14-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
% -19  2-15-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
% -20  2-16-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat 15 x 1 x 6000/10000 ms B
APcrit = [-7 -1 1 -6]; Nav = 1; % decent thr def
for ifreq=1:15,
    JbB_70dB(ifreq) = JLbeat('RG10216b', '2-1-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
    JbC_70dB(ifreq) = JLbeat('RG10216b', '2-2-BFS', ifreq, Nav, APcrit, 0, 'C_70dB');
    JbI_70dB(ifreq) = JLbeat('RG10216b', '2-3-BFS', ifreq, Nav, APcrit, 0, 'I_70dB');
    %
    JbB_60dB(ifreq) = JLbeat('RG10216b', '2-4-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10216b', '2-5-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10216b', '2-6-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB2_60dB(ifreq) = JLbeat('RG10216b', '2-13-BFS', ifreq, Nav, APcrit, 0,'B2_60dB');
    JbB3_60dB(ifreq) = JLbeat('RG10216b', '2-14-BFS', ifreq, Nav, APcrit, 0, 'B3_60dB');
    JbB4_60dB(ifreq) = JLbeat('RG10216b', '2-15-BFS', ifreq, Nav, APcrit, 0, 'B4_60dB');
    JbB5_60dB(ifreq) = JLbeat('RG10216b', '2-16-BFS', ifreq, Nav, APcrit, 0, 'B5_60dB');
    %
    JbB_50dB(ifreq) = JLbeat('RG10216b', '2-7-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10216b', '2-8-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10216b', '2-9-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    %
    JbB_40dB(ifreq) = JLbeat('RG10216b', '2-10-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbI_40dB(ifreq) = JLbeat('RG10216b', '2-11-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
    JbC_40dB(ifreq) = JLbeat('RG10216b', '2-12-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
end
local_addComments('mini spikes, large EPSPS', 'decent thr def');
save(SFN, 'Jb*'); clear Jb*;


function local_RG10219_1(SFN);
% RG10219 cell 1 ("spikes & EPSPS")  
%   cell 1 ("spikes & EPSPs")
%      30 dB contra/ipsi
%      40 dB bin/ipsi/contra
%      50 dB bin/ipsi/contra
%      60 dB bin
%      70 dB bin
% -1   1-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat               15 x 1 x 6000/10000 ms B
% -2   1-2-BFS  200*:50*:800 Hz   70 dB     4 Hz beat               13 x 1 x 6000/10000 ms B
% -3   1-3-BFS  200*:50*:800 Hz   50 dB     4 Hz beat               13 x 1 x 6000/10000 ms B
% -4   1-4-BFS  200*:50*:800 Hz   50|-50 dB 4 Hz beat               13 x 1 x 6000/10000 ms B
% -5   1-5-BFS  200*:50*:800 Hz   -50|50 dB 4 Hz beat               13 x 1 x 6000/10000 ms B
% -6   1-6-BFS  200*:50*:800 Hz   -40|40 dB 4 Hz beat               13 x 1 x 6000/10000 ms B
% -7   1-7-BFS  200*:50*:800 Hz   40 dB     4 Hz beat               13 x 1 x 6000/10000 ms B
% -8   1-8-BFS  200*:50*:800 Hz   40|-40 dB 4 Hz beat               13 x 1 x 6000/10000 ms B
% -9   1-9-BFS  200*:50*:800 Hz   30|-30 dB 4 Hz beat               13 x 1 x 6000/10000 ms B
% -10  1-10-BFS 200*:50*:800 Hz   -30|30 dB 4 Hz beat               13 x 1 x 6000/10000 ms B
APcrit = [-8 -1.5 1.5 -6];  Nav = 1; % decent thr def
for ifreq=1:15, % 100:100:1500 Hz
    JbB_60dB(ifreq) = JLbeat('RG10219', '1-1-BFS', ifreq,3, APcrit, 0, 'B_60dB'); 
end
for ifreq=1:13, % 200:50:800 Hz
    JbB_70dB(ifreq) = JLbeat('RG10219', '1-2-BFS', ifreq, Nav, APcrit, 0, 'B_70dB'); 
    JbB_50dB(ifreq) = JLbeat('RG10219', '1-3-BFS', ifreq, Nav, APcrit, 0, 'B_50dB'); 
    JbI_50dB(ifreq) = JLbeat('RG10219', '1-4-BFS', ifreq, Nav, APcrit, 0, 'I_50dB'); 
    JbC_50dB(ifreq) = JLbeat('RG10219', '1-5-BFS', ifreq, Nav, APcrit, 0, 'C_50dB'); 
    JbC_40dB(ifreq) = JLbeat('RG10219', '1-6-BFS', ifreq, Nav, APcrit, 0, 'C_40dB'); 
    JbB_40dB(ifreq) = JLbeat('RG10219', '1-7-BFS', ifreq, Nav, APcrit, 0, 'B_40dB'); 
    JbI_40dB(ifreq) = JLbeat('RG10219', '1-8-BFS', ifreq, Nav, APcrit, 0, 'I_40dB'); 
    JbI_30dB(ifreq) = JLbeat('RG10219', '1-9-BFS', ifreq, Nav, APcrit, 0, 'I_30dB'); 
    JbC_30dB(ifreq) = JLbeat('RG10219', '1-10-BFS', ifreq, Nav, APcrit, 0, 'C_40dB'); 
end
local_addComments('spikes & EPSPS', 'decent thr def');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10219_2(SFN);
% RG10219 cell 2 ("mini spikes, large EPSPS")  
%      40 dB bin
%      50 dB bin/ipsi/contra
%      60 dB bin
% -11  2-1-BFS  100*:100*:1400 Hz 50 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -12  2-2-BFS  100*:100*:1400 Hz 60 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -13  2-3-BFS  100*:100*:1400 Hz 40 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -14  2-4-BFS  100*:100*:1400 Hz 50|-50 dB 4 Hz beat               14 x 1 x 6000/10000 ms B
% -15  2-5-BFS  100*:100*:1400 Hz -50|50 dB 4 Hz beat               14 x 1 x 6000/10000 ms B
APcrit = [-5 -1 1 -4]; Nav = 1;
for ifreq=1:14,
    JbB_50dB(ifreq) = JLbeat('RG10219', '2-1-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbB_60dB(ifreq) = JLbeat('RG10219', '2-2-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbB_40dB(ifreq) = JLbeat('RG10219', '2-3-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbI_50dB(ifreq) = JLbeat('RG10219', '2-4-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10219', '2-5-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
end
local_addComments('mini spikes, large EPSPS', 'poorly def thr; hardly spiking, but binaural waveforms');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10219_3(SFN);
% RG10219 cell 3 ("spikes + EPSPS", weak binaurality; mostly ipsi driven; very few spikes. -10 V/s is correct thr)  
% JL XLS: cell was before turnover range, but was little bit phase sensitive
%      30 dB bin
%      40 dB bin
%      50 dB bin
%      60 dB bin/ipsi/contra
%      70 dB bin
%      80 dB bin/ipsi/contra
% -16  3-1-BFS  100*:100*:1400 Hz 50 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -17  3-2-BFS  100*:100*:1400 Hz 40 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -18  3-3-BFS  100*:100*:1400 Hz 70 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -19  3-4-BFS  100*:100*:1400 Hz 60 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -20  3-5-BFS  100*:100*:1400 Hz 60|-60 dB 4 Hz beat               14 x 1 x 6000/10000 ms B
% -21  3-6-BFS  100*:100*:1400 Hz -60|60 dB 4 Hz beat               14 x 1 x 6000/10000 ms B
% -41  3-26-BFS 100*:100*:1400 Hz 80 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -42  3-27-BFS 100*:100*:1400 Hz 80|-80 dB 4 Hz beat               14 x 1 x 6000/10000 ms B
% -43  3-28-BFS 100*:100*:1400 Hz -80|80 dB 4 Hz beat               14 x 1 x 6000/10000 ms B
% -44  3-29-BFS 100*:100*:1400 Hz 30 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
APcrit = [-10 -1 1.4]; Nav = 1;
for ifreq=1:14,
    JbB_50dB(ifreq) = JLbeat('RG10219', '3-1-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbB_40dB(ifreq) = JLbeat('RG10219', '3-2-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbB_70dB(ifreq) = JLbeat('RG10219', '3-3-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
    JbB_60dB(ifreq) = JLbeat('RG10219', '3-4-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10219', '3-5-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
    JbC_60dB(ifreq) = JLbeat('RG10219', '3-6-BFS', ifreq, Nav, APcrit, 0, 'C_60dB');
    JbB_80dB(ifreq) = JLbeat('RG10219', '3-26-BFS', ifreq, Nav, APcrit, 0, 'B_80dB');
    JbI_80dB(ifreq) = JLbeat('RG10219', '3-27-BFS', ifreq, Nav, APcrit, 0, 'I_80dB');
    JbC_80dB(ifreq) = JLbeat('RG10219', '3-28-BFS', ifreq, Nav, APcrit, 0, 'C_80dB');
    JbB_30dB(ifreq) = JLbeat('RG10219', '3-29-BFS', ifreq, Nav, APcrit, 0, 'B_30dB');
end
local_addComments('spikes + EPSPS', 'weak binaurality; mostly ipsi driven; very few spikes. well def thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10219_4(SFN);
% RG10219 cell 4 ("mini spikes and big EPSPs" )  
%      40 dB bin
%      50 dB bin
%      60 dB bin
%      70 dB bin
% -45  4-1-BFS  100*:100*:1400 Hz 60 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -46  4-2-BFS  100*:100*:1400 Hz 50 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -47  4-3-BFS  100*:100*:1400 Hz 40 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -48  4-4-BFS  100*:100*:1400 Hz 70 dB     4 Hz beat               14 x 1 x 6000/10000 ms B

APcrit = [-9 -1 1.4 -7]; Nav = 1; % reasonable thr
for ifreq=1:14,
    JbB_60dB(ifreq) = JLbeat('RG10219', '4-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbB_50dB(ifreq) = JLbeat('RG10219', '4-2-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbB_40dB(ifreq) = JLbeat('RG10219', '4-3-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbB_70dB(ifreq) = JLbeat('RG10219', '4-4-BFS', ifreq, Nav, APcrit, 0, 'B_70dB');
end
local_addComments('mini spikes and big EPSPs', 'reasonable thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10219_5(SFN);
% RG10219 cell 5 ("mini spikes and big EPSPs" )  
%      30 dB bin
%      40 dB bin
%      50 dB bin
%      60 dB bin
% -49  5-1-BFS  100*:100*:1400 Hz 60 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -50  5-2-BFS  100*:100*:1400 Hz 40 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -51  5-3-BFS  100*:100*:1400 Hz 50 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
% -52  5-4-BFS  100*:100*:1400 Hz 30 dB     4 Hz beat               14 x 1 x 6000/10000 ms B
APcrit = [-6 -1 1.4]; Nav = 1; 
for ifreq=1:14,
    JbB_60dB(ifreq) = JLbeat('RG10219', '5-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbB_40dB(ifreq) = JLbeat('RG10219', '5-2-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbB_50dB(ifreq) = JLbeat('RG10219', '5-3-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbB_30dB(ifreq) = JLbeat('RG10219', '5-4-BFS', ifreq, Nav, APcrit, 0, 'B_30dB');
end
local_addComments('mini spikes and big EPSPs', 'reasonable thr; good S/N');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10219_6(SFN);
% RG10219 cell 6 ("spikes and EPSPs" )  
%      60 dB bin/ipsi
% -53  6-1-BFS  100*:100*:1200 Hz 60 dB     4 Hz beat               12 x 1 x 6000/10000 ms B
% -54  6-2-BFS  100*:100*:1200 Hz 60|-60 dB 4 Hz beat               12 x 1 x 6000/10000 ms B
APcrit = [-8 -1 1.4 -6]; Nav = 1; % okay thr
for ifreq=1:12,
    JbB_60dB(ifreq) = JLbeat('RG10219', '6-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbI_60dB(ifreq) = JLbeat('RG10219', '6-2-BFS', ifreq, Nav, APcrit, 0, 'I_60dB');
end
local_addComments('spikes and EPSPs', 'okay thr');
save(SFN, 'Jb*'); clear Jb*;


function local_RG10223_1(SFN);
% RG10223 cell 1 ("spikes and EPSPs", "Gain was only 10, changed to 50 after this recording" )  
%      60 dB bin
% -3   1-3-BFS 100*:100*:2500 Hz 60 dB     4 Hz beat 25 x 1 x 6000/9500 ms B
APcrit = [-8 -1 1.4]; Nav = 1; % need high thr to correct quantization noise
for ifreq=1:25,
    JbB_60dB(ifreq) = JLbeat('RG10223', '1-3-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
end
local_addComments('spikes and EPSPs; Gain was only 10, changed to 50 after this recording', 'need high thr to correct quantization noise');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10223_2(SFN);
% RG10223 cell 2 ("spikes + EPSPs". "isolation gets worse (2-2 and up)")  
%      50 dB 2xbin/ipsi/contra
% -4   2-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -5   2-2-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -6   2-3-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -7   2-4-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
APcrit = [-9 -1 1.4 -7]; Nav = 1; % okay thr
for ifreq=1:15,
    JbB_50dB(ifreq) = JLbeat('RG10223', '2-1-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10223', '2-2-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10223', '2-3-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbB2_50dB(ifreq) = JLbeat('RG10223', '2-4-BFS', ifreq, Nav, APcrit, 0, 'B2_50dB');
end
local_addComments('spikes + EPSPs; isolation gets worse (2-2 and up)', 'okay thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10223_3(SFN);
% RG10223 cell 3 ("spikes + EPSPs")  
%      50 dB 2xbin/ipsi/contra
% -8   3-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -9   3-2-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -10  3-3-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -11  3-4-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
APcrit = [-7 -1 1.4 -6]; Nav = 1; % okay thr
for ifreq=1:15,
    JbB_50dB(ifreq) = JLbeat('RG10223', '3-1-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10223', '3-2-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10223', '3-3-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbB2_50dB(ifreq) = JLbeat('RG10223', '3-4-BFS', ifreq, Nav, APcrit, 0, 'B2_50dB');
end
local_addComments('spikes + EPSPs', 'okay thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10223_4(SFN);
% RG10223 cell 4 ("spikes + EPSPs")  
%      50 dB 2xbin/ipsi/contra
%      40 dB bin/ipsi
% -12  4-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -13  4-2-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -14  4-3-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -15  4-4-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -16  4-5-BFS 100*:100*:1500 Hz 40 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -17  4-6-BFS 100*:100*:1500 Hz 40|-40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
APcrit = [-7 -1 1.4 -6]; Nav = 1; % not a very crisp thr
for ifreq=1:15,
    JbB_50dB(ifreq) = JLbeat('RG10223', '4-1-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10223', '4-2-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10223', '4-3-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbB2_50dB(ifreq) = JLbeat('RG10223', '4-4-BFS', ifreq, Nav, APcrit, 0, 'B2_50dB');
    JbB_40dB(ifreq) = JLbeat('RG10223', '4-5-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbI_40dB(ifreq) = JLbeat('RG10223', '4-6-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
end
local_addComments('spikes + EPSPs', 'not a very crisp thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10223_5(SFN);
% RG10223 cell 5 ("spikes + EPSPs, but losing cells slowly after 1000 Hz", isolation getting worse)  
%      50 dB bin
%      40 dB bin/ipsi/contra
% -18  5-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -19  5-2-BFS 100*:100*:1500 Hz 40 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -20  5-3-BFS 100*:100*:1500 Hz -40|40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -21  5-4-BFS 100*:100*:1500 Hz 40|-40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
APcrit = [-6.5 -1 1.4 -5.5]; Nav = 1; % not a very crisp separation
for ifreq=1:15,
    JbB_50dB(ifreq) = JLbeat('RG10223', '5-1-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbB_40dB(ifreq) = JLbeat('RG10223', '5-2-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbC_40dB(ifreq) = JLbeat('RG10223', '5-3-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
    JbI_40dB(ifreq) = JLbeat('RG10223', '5-4-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
end
local_addComments('spikes + EPSPs, but losing cells slowly after 1000 Hz; isolation getting worse', 'not a very crisp thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10223_6(SFN);
% RG10223 cell 6 ("recordings after a Gigaseal, R=40Mohm, not standard waveforms")  
%      20 dB bin
%      30 dB bin
%      40 dB bin/ipsi/contra
%      50 dB bin
%      60 dB bin
%      80 dB bin
% -22  6-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -23  6-2-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -24  6-3-BFS 100*:100*:1500 Hz 40 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -25  6-4-BFS 100*:100*:1500 Hz 40|-40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -26  6-5-BFS 100*:100*:1500 Hz -40|40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -27  6-6-BFS 100*:100*:1500 Hz 30 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -28  6-7-BFS 100*:100*:1500 Hz 20 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% XXX -29  6-8-BFS 100*:100*:1500 Hz 80 dB     4 Hz beat 1* x 1 x 6000/9500 ms B aborted due to Voltage clamp
% -30  6-9-BFS 100*:100*:1500 Hz 80 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
APcrit = [-4.5 -1 1.4 -3.75]; Nav = 1; % not a very crisp thr
for ifreq=1:15,
    JbB_60dB(ifreq) = JLbeat('RG10223', '6-1-BFS', ifreq, Nav, APcrit, 0, 'B_60dB');
    JbB_50dB(ifreq) = JLbeat('RG10223', '6-2-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbB_40dB(ifreq) = JLbeat('RG10223', '6-3-BFS', ifreq, Nav, APcrit, 0, 'B_40dB');
    JbI_40dB(ifreq) = JLbeat('RG10223', '6-4-BFS', ifreq, Nav, APcrit, 0, 'I_40dB');
    JbC_40dB(ifreq) = JLbeat('RG10223', '6-5-BFS', ifreq, Nav, APcrit, 0, 'C_40dB');
    JbB_30dB(ifreq) = JLbeat('RG10223', '6-6-BFS', ifreq, Nav, APcrit, 0, 'B_30dB');
    JbB_20dB(ifreq) = JLbeat('RG10223', '6-7-BFS', ifreq, Nav, APcrit, 0, 'B_20dB');
    JbB_80dB(ifreq) = JLbeat('RG10223', '6-9-BFS', ifreq, Nav, APcrit, 0, 'B_80dB');
end
local_addComments('recordings after a Gigaseal, R=40Mohm, not standard waveforms', 'not a very crisp thr');
save(SFN, 'Jb*'); clear Jb*;

function local_RG10223_7(SFN);
% RG10223 cell 7 ("failed patch attempt: small upward spikes", "lots of inhibition, not clear beats")  
%      50 dB 2xbin/ipsi/contra
% -31  7-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -32  7-2-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -33  7-3-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -34  7-4-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
APcrit = [-5 -1 1.4]; Nav = 1; % not a well defined thr
for ifreq=1:15,
    JbB_50dB(ifreq) = JLbeat('RG10223', '7-1-BFS', ifreq, Nav, APcrit, 0, 'B_50dB');
    JbB2_50dB(ifreq) = JLbeat('RG10223','7-2-BFS', ifreq, Nav, APcrit, 0, 'B2_50dB');
    JbI_50dB(ifreq) = JLbeat('RG10223', '7-3-BFS', ifreq, Nav, APcrit, 0, 'I_50dB');
    JbC_50dB(ifreq) = JLbeat('RG10223', '7-4-BFS', ifreq, Nav, APcrit, 0, 'C_50dB');
end
local_addComments('failed patch attempt: small upward spikes', 'not a well defined thr');
save(SFN, 'Jb*'); clear Jb*;


function  SL = local_slist(S);
% cnvert struc to struct array
FNS = fieldnames(S);
SL = [];
for ii=1:numel(FNS),
    fn = FNS{ii};
    iCell = S.(fn);
    for jj=1:numel(iCell),
        SL = [SL, struct('ExpID', fn, 'icell', iCell(jj))];
    end
end
 






