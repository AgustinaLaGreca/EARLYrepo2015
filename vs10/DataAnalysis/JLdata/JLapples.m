function JLapples
% JL apples - 2nd round of JL data processing: JLbeatStats, PCA, etc

% todo:
%    - HP filtering spont activity (STD)
%    - AP truncation & HP filtering stim act (STD)
%    D bookkeeping: link monaural data to binaural data -> linearity analysis
%    - estimate # inputs
%    - automated junk rejection
%    - 
%    - compare linear pred to real-time recording
%    - database: peak counts of spont act; use peak sharpness.
%    - compare peak dynamics and try to understand in terms of
%      synchronization across inputs.
%    - monaural data: compare baseline variance to stimulated variance. Is
%      "spontaneous variance" added to entrained activity? Compare to 
%      binaural case, where both inputs are phase locked. Lookat short term
%      variance.
%    - aim: analyze nonlinear binaural interaction. Escape the
%      linearization induced by periodic averaging.


% Check AP thr of RG10191-15 <3-3-BFS> 100 Hz 

load D:\processed_data\JL\JLbeatStats\STall.mat; % This loads ST previous returned by JLbeatstats('all')

%=======multi peaked stuff==============
iok = ([ST.Freq1]<=300) & ([ST.AC]>0.5); % low stim freq & decent S/N ratio
STm = ST(iok);
iok = ([ST.Freq1]<=300) & ([ST.StimType]~='B') & ([ST.AC]>0.5); % low stim freq & monaural stim & decent S/N ratio
STmon = ST(iok);

P = JLprincomp(STm(1),1); % RG10191 B2_60dB(3):  200 Hz, [60 60] dB
P = JLprincomp(STm(2),1); % RG10191 B2_60dB(4):  250 Hz, [60 60] dB
P = JLprincomp(STm(3),1); % RG10191 B2_60dB(5):  300 Hz, [60 60] dB
P = JLprincomp(STm(4),1); % RG10191 B_60dB(3):  200 Hz, [60 60] dB
P = JLprincomp(STm(5),1); % RG10191 B_60dB(4):  250 Hz, [60 60] dB
P = JLprincomp(STm(6),1); % RG10191 B_60dB(5):  300 Hz, [60 60] dB
P = JLprincomp(STm(7),1); % RG10191 C_60dB(3):  200 Hz, [-60 60] dB
P = JLprincomp(STm(8),1); % RG10191 C_60dB(4):  250 Hz, [-60 60] dB
P = JLprincomp(STm(9),1); % RG10191 C_60dB(5):  300 Hz, [-60 60] dB

P = JLprincomp(STm(19),1); % interesting Ipsi

P = JLprincomp(STm(23),1); % RG10197 B2_70dB(3):  200 Hz, [70 70] dB
P = JLprincomp(STm(24),1); % RG10197 B2_70dB(5):  300 Hz, [70 70] dB

% =========mon===============
P = JLprincomp(STmon(1),1); 

P = JLprincomp(STmon(4),1); % two inputs?

P = JLprincomp(STmon(6),1);

% searching in vain for indep inputs
iii=0;P = JLprincomp(STmon(9),1); c = P.C2; ibig = find(c.pCF(:,1)<-0.0171); ismall = find(c.pCF(:,1)>-0.0171);
% RG10198 C_70dB(2):  100 Hz, [-70 70] dB
P = JLprincomp(STmon(9),1, {[] ibig});
P = JLprincomp(STmon(9),1, {[] ismall});


%%%==========
% two components??
P = JLprincomp(STmon(13),1, {[] []}); % RG10216b C_70dB(2):  200 Hz, [-70 70] dB


% ==========search for significant multpeak data=================
% contra stim
ihitc = find([ST.Vstim_STD]./[ST.Vspont_STD]>1.2 & [ST.StimType]=='C' & [ST.AC]>0.25 & [ST.Freq1]<400);
STc = ST(ihitc);
hist([STc.Vcontra_Npeak], 1:max([STc.Vcontra_Npeak]))
find([STc.Vcontra_Npeak]>2)
%   1     2     3     6     8     9    15    17    18    19    21    23    24    25    26    27    30    31    34    36
%     37    38    41    42    43    45    47    48    49    53    54    55    56
find([STc.Vcontra_Npeak]>3)
%    9    17    18    21    36    42    43    47    48    49    53
%
% ipsi stim
ihiti = find([ST.Vstim_STD]./[ST.Vspont_STD]>1.2 & [ST.StimType]=='I' & [ST.AC]>0.25 & [ST.Freq1]<400);
STi = ST(ihiti);
hist([STi.Vipsi_Npeak], 1:max([STi.Vipsi_Npeak]))
find([STi.Vipsi_Npeak]>2)
%   3     4     5     6     8     9    13    14    15    16    17    21    23

[P iselect]=JLprincomp(STc(17),1,[],[0.2 1])
JLprincomp(STc(17),1,[],[0 0.2])

%============== clickable scatter plot of waveform STDs=================
hl=plot([s.Vspont_STD], [s.Vstim_STD], '.');
IDpoints(hl,s ,1:numel(s), @(s,ii)[num2str(ii) ':  ' s(ii).TTT], ...
    'raw data', @(s,ii)clampfit(s(ii)), ...
    'periodicity analysis', @(s,ii)JLbeatPlot(s(ii)), ...
    'display stats', @(s,ii)disp(s(ii)), ...
    'PCA', @(s,ii)JLprincomp(s(ii)), ...
    'SpontPeaks', @(s,ii)JLspontPeaks(s(ii)) ...
    );


%---------spont activity studies per cell-----------
allCellIdx = unique([ST.UniqueCellIndex]);
ii = 0;
% repeat next line
ii=ii+1, STcl = ST([ST.UniqueCellIndex]==allCellIdx(ii)); JLspontAna(STcl)

%     ExpID: 'RG10195' cell 3
STcl = ST([ST.UniqueCellIndex]==allCellIdx(4)); JLspontAna(STcl); % very oscillatory
JLreadBeats(STcl(1).ExpID, STcl(1).icell)
% very multipeaked:
% 100 Hz  R = [0 0 0] 0 spikes; rho=0.28;
% 1: 100 Hz  RG10195-48 <3-10-BFS> 0.1*:0.1*:1 kHz -65|65 dB 4 Hz beat 10 x 1 x 6000/8500 ms B  <c164.abf> 
JLprincomp(JbC_65dB(1));

P = JLprincomp(JbC_75dB(1), 1, [], 0.8);
f4; xdplot(P.C2.dt, mean(P.C2.rec,2), 'k', 'linewidth', 8)

% pretty oscillatory
STcl = ST([ST.UniqueCellIndex]==allCellIdx(5)); JLspontAna(STcl);
JLreadBeats(STcl(1).ExpID, STcl(1).icell)
JLprincomp(JbC_80dB(2),1,[],0.7); % what's this all-or-nothing behavior? Peaks are 4 ms apart!


%===========================================
% beat variance analysis
STipsi_100Hz = ST([ST.StimType]=='I' & [ST.Freq1]==100);
STcontra_100Hz = ST([ST.StimType]=='C' & [ST.Freq1]==100);

JLbeatVar(STipsi_100Hz(7)); % in between spikes, the variance is below spont value

JLbeatVar(ST(176)); % 1: 100 Hz  RG10195-48 <3-10-BFS> 0.1*:0.1*:1 kHz -65|65 dB
JLprincomp(ST(176));

























