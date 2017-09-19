function D = ZWOAEphase(FN, varargin);
% ZWOAEphase - phase analysis of ZWOA data
%   ZWOAEphase('foo') analyzes phase data in file foo.
%
%   ZWOAEphase('foo', 15) analyzes phase data in file foo015.
%
%
%  See also getZWOAEdata.

% read data
D = getZWOAEdata(FN, varargin{:});
% extract important stimulus params
sp = D.stimpars;
i1 = [sp.N1~=0] & [sp.L1>0] ;
i2 = [sp.N2~=0] & [sp.L2>0] ;
F1 = sp.F1(i1).'; % lower freq(s) in Hz
F2 = sp.F2(i2).'; % upper freq(s) in Hz
dt = 1/D.fs; % sample period in ms
Nsam = length(D.signal);
FmaxPlot = min(2.2*max(F1), 15e3);

% by convention F1 is the zwuis complex, even if it's the upper freq
if length(F2)>1, [F1, F2]=swap(F1, F2); end
Nzw  = length(F1); % # tones in F1
Fprim = [F1; F2]; % col vector containing all the primaries

PS = ZWOAEplotStruct; % universal worldwide ZWOAEzwiep plotting conventions

% spectrum
df = 1e3/(Nsam*dt); % spectral spacing in Hz
freq = df*(0:Nsam-1);
Cspec = fft(D.signal).'; % complex spectrum in row vector
PH = angle(Cspec)/2/pi; % phase in cycles
MG = A2dB(abs(Cspec)*2/Nsam); % magnitude in dB
MG = MG-D.micgain-D.micsensitivity.dBV+D.micsensitivity.dBSPL;

% Primaries
Mprim = local_SpecComp(df, MG, Fprim);
Pprim = local_SpecComp(df, PH, Fprim);

% plot magn spectrum
set(gcf,'position', [12 47 565 687])
subplot(2,1,1);
plot(freq,MG, PS.Spec);
xplot(freq(2:2:end),MG(2:2:end),PS.Floor); % noise floor: uneven components containing no stimulus periodicity
xplot(Fprim, Mprim, PS.Prim);
xlim([0 FmaxPlot]); ylim(max(MG) +[-90 5]);
title([FN, ' ' num2str(D.index) ': #speakers=' num2str(D.setupinfo.numchan)],'interpreter', 'none')

% basic matrices
DM = EvaDifMat(Nzw); % matrix for evaluating differences across F1 cmps
SM = abs(DM); % matrix for evaluating sums of unequal F1 cmps
SMp = [SM; 2*eye(Nzw)]; % extend SM to include equal pairs

LG = {}; % legend cellstring

% =========F1+F1-F2 terms========
Nsum = size(SMp,1);
Mc1 = [SMp, -ones(Nsum,1)]; % the -ones always subtract the last prim cmp, ie F2
[Fc1, MGc1, PHc1, Rdef]=local_DP(df, MG, PH, Fprim, Mc1); Rdef
measMGc1 = local_SpecComp(df, MG, Fc1); % measured levels of DPs
subplot(2,1,1);
xplot(Fc1, measMGc1,PS.DPnear);
%xplot(Fc1, mean(measMGc1)+MGc1-mean(MGc1),'k.');
subplot(2,1,2);
plot(Fc1, PHc1,PS.DPnear); xlim([0 FmaxPlot]); grid on
LG = [LG local_delayStr(Fc1, PHc1)];
Mc1_r = [ones(1,length(Fprim)); Mc1]; % provide extra row of ones to lift rank deficiency
PH1c1 = Mc1_r\[0; PHc1];
xplot(Fprim,PH1c1,PS.DPnear);
LG = [LG local_delayStr(Fprim, PH1c1,1)];
PHRecalc_c1 = Mc1*PH1c1; 
xplot(Fc1, PHRecalc_c1, '.k')
LG = [LG local_delayStr(Fc1, PHRecalc_c1)];

% =========lower F1-F1+F2 terms========
Ndif = size(DM,1);
Mslo = [DM, ones(Ndif,1)]; % the ones always add the last prim cmp, ie F2
[Fslo, MGslo, PHslo, Rdef]=local_DP(df, MG, PH, Fprim, Mslo); Rdef
measMGslo = local_SpecComp(df, MG, Fslo); % measured levels of DPs
subplot(2,1,1);
xplot(Fslo, measMGslo,PS.DPsuplo);
%xplot(Fslo, mean(measMGslo)+MGslo-mean(MGslo),'k.');
subplot(2,1,2);
xplot(Fslo, PHslo,PS.DPsuplo); 
LG = [LG local_delayStr(Fslo, PHslo)];
Mslo_r = [ones(1,length(Fprim)); Mslo]; % provide extra row of ones to lift rank deficiency
PH1slo = Mslo_r\[0; PHslo];
xplot(Fprim,PH1slo,PS.DPsuplo);
LG = [LG local_delayStr(Fprim,PH1slo,1)];
PHRecalc_slo = Mslo*PH1slo; 
xplot(Fslo, PHRecalc_slo, '.k')
LG = [LG local_delayStr(Fslo, PHRecalc_slo)];

% =========upper F1-F1+F2 terms========
Mshi = [-DM, ones(Ndif,1)]; % the ones always add the last prim cmp, ie F2
[Fshi, MGshi, PHshi, Rdef]=local_DP(df, MG, PH, Fprim, Mshi); Rdef
measMGshi = local_SpecComp(df, MG, Fshi); % measured levels of DPs
subplot(2,1,1);
xplot(Fshi, measMGshi,PS.DPsuphi);
%xplot(Fshi, mean(measMGshi)+MGshi-mean(MGshi),'k.');
subplot(2,1,2);
xplot(Fshi, PHshi,PS.DPsuphi); 
LG = [LG local_delayStr(Fshi, PHshi)];
Mshi_r = [ones(1,length(Fprim)); Mshi]; % provide extra row of ones to lift rank deficiency
PH1shi = Mshi_r\[0; PHshi];
xplot(Fprim,PH1shi,PS.DPsuphi);
LG = [LG local_delayStr(Fprim,PH1shi,1)];
PHRecalc_shi = Mshi*PH1shi; 
xplot(Fshi, PHRecalc_shi, '.k')
LG = [LG local_delayStr(Fshi, PHRecalc_shi)];

% =========ALL F1-F1+F2 terms========
Msup = [[DM; -DM], ones(2*Ndif,1)]; % the ones always add the last prim cmp, ie F2
[Fsup, MGsup, PHsup, Rdef]=local_DP(df, MG, PH, Fprim, Msup); Rdef
measMGsup = local_SpecComp(df, MG, Fsup); % measured levels of DPs
subplot(2,1,1);
%xplot(Fsup, measMGsup,'co', 'markerfacecolor', 'c', 'markersize', 4);
%xplot(Fsup, mean(measMGsup)+MGsup-mean(MGsup),'.', 'color', 0.7*[1 1 1]);
subplot(2,1,2);
%xplot(Fsup, PHsup,'c*'); 
Msup_r = [ones(1,length(Fprim)); Msup]; % provide extra row of ones to lift rank deficiency
PH1sup = Msup_r\[0; PHsup];
xplot(Fprim,PH1sup,PS.DPsupall);
LG = [LG local_delayStr(Fprim,PH1sup,1)];
PHRecalc_sup = Msup*PH1sup; 
xplot(Fsup, PHRecalc_sup, '.', 'color', 0.7*[1 1 1])
LG = [LG local_delayStr(Fsup, PHRecalc_sup)];

% =========2*F2-F1 terms========
Mc2 = [-eye(Nzw), 2*ones(Nzw,1)];
[Fc2, MGc2, PHc2, Rdef]=local_DP(df, MG, PH, Fprim, Mc2); Rdef
MGc2 = local_SpecComp(df, MG, Fc2); % measured levels of DPs
subplot(2,1,1);
xplot(Fc2, MGc2, PS.DPfar);
subplot(2,1,2);
xplot(Fc2, PHc2, PS.DPfar);
LG = [LG local_delayStr(Fc2, PHc2)];

subplot(2,1,2);
legend(LG{:}, 'location', 'West');




%==========LOCALS==================================
function [Fdp, Mdp, Pdp, Rdef]=local_DP(df, MG, PH, Fprim, DPmat);
% Compute DP frequencies and DP reference phases & magnitudes from primaries
% using an interaction matrix.
% Input args:
%      df: frequency spacing in Hz
%      MG: magnitude spectrum (running from 0 kHz) in dB
%      PH: phase spectrum (running from 0 kHz) in Cycles
%   Fprim: primary frequencies (in Hz) in column vector
%   DPmat: interaction matrix; DPmat*Fprim yield DP frequencies
% Output args:
%     Fdp: (unsorted) frequencies (in Hz) of distortions
%     Mdp: magn of distortions in dB expressed re primary magnitudes
%     Pdp: phase of distortions in cycles expressed re primary phases
Mprim = local_SpecComp(df, MG, Fprim); % Primary magnitudes
Pprim = local_SpecComp(df, PH, Fprim); % Primary phases
Fdp = DPmat*Fprim; % DP freqs
Mref = DPmat*Mprim; % reference magnitudes @ dp frequencies as derived from primary magnitudes
Pref = DPmat*Pprim; % reference phases @ dp frequencies as derived from primary phases
Mdp = local_SpecComp(df, MG, Fdp)-Mref; % dp amgnitudes expressed re primary magnitues
Pdp = local_SpecComp(df, PH, Fdp)-Pref; % idem phases
% unwrap phases along freq axis. Note: cannot sort the DPs since they are
% needed for leftdivision using DPmat when estimating underlying primary
% contributions.
[dum, isort] = sort(Fdp);
Pdp = cunwrap(Pdp(isort));
Pdp = Pdp - round(mean(Pdp)); % move to the cycle closest to zero
Pdp(isort) = Pdp; % unsort
Rdef = size(DPmat,2)-rank(DPmat); % rank deficiency of the matrix

function S = local_SpecComp(df, Spec, Freq);
% return components of spectrum S @ requested frequencies;
ii = 1+round(Freq/df); % indices of components in spectrum
S = Spec(ii);

function Str = local_delayStr(fr, ph, ignoreLast);
% estimate group delay from freq/phase data
if nargin<3, ignoreLast = 0; end
if ignoreLast,
    [fr, ph] = deal(fr(1:end-1), ph(1:end-1));
end
qq = polyfit(fr,ph,1);
Str = [num2str(round(-1e6*qq(1))) ' \mus'];







