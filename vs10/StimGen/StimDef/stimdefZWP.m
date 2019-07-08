function Params = stimdefZWP(EXP);
% stimdefZWP - definition of stimulus and GUI for ZWP stimulus paradigm
%    P=stimdefZWP(EXP) returns the definition for the ZWP (zwuis with profile)
%    ZWP is di
%
%    See also stimGUI, stimDefDir, Experiment, makestimZW.

PairStr = ' Pairs of numbers are interpreted as [left right].';
ClickStr = ' Click button to select ';

ZW = local_zwuisPanel;
PR = local_profilePanel;
SPL = SPLpanel('-', EXP, '' , 'Component');
%Dur = DurPanel('-', EXP);
Pres = local_presPanel;
%====================
Params=GUIpiece('Params'); % upper half of GUI: parameters
Params = add(Params, ZW, 'below', [20 0]);
Params = add(Params, PR, below(ZW));
Params = add(Params, SPL, nextto(ZW), [20 0]);
Params = add(Params, Pres, below(SPL) ,[30 10]);
Params = add(Params, PlayTime(), below(Pres) ,[-75 7]);

%=================================
function P = local_zwuisPanel;
LowFreq = ParamQuery('LowFreq', 'low:', '15000.5', 'Hz', ...
    'rreal/positive', 'Approximate lowest frequency of zwuis complex.', 1);
HighFreq = ParamQuery('HighFreq', 'high:', '15000.5', 'Hz', ...
    'rreal/positive', 'Approximate highest frequency of zwuis complex.', 1);
Ncomp = ParamQuery('Ncomp', '#comp', '123', '', ...
    'posint', 'Number of tones in complex.' ,1);
BaseFreq = ParamQuery('BaseFreq', 'base:', '0.128 ', 'Hz', ...
    'rreal/positive', 'Base frequency. Every component frequency will be an integer multiple of this freq.', 1);
FreqTol = ParamQuery('FreqTol', 'maxdev:', '120', '%', ...
    'rreal/positive', ['Maximum deviation of component frequencies from regular spacing, ' char(10), ...
    'expressed in % of average component spacing.'], 1);
% Rseed = paramquery('Rseed', 'seed:', '8445963002', '', ...
%     'rseed', 'Random seed used for frequencies and phases. Specify NaN to refresh seed upon each realization.',1);
Nit = ParamQuery('Nit', '#iter:', '123', '', ...
    'posint', 'Max number of iterations to realize zwuis frequencies.' ,1);
ZWseed = ParamQuery('ZWseed', 'seed:', '12345678', '', ...
    'rseed', 'Random seed used to realize zwuis frequencies. Nan means refresh every time, EXP or EXP3, etc, means one seed for each experiment.' ,1);
ZwuisMess = messenger('ZwuisMess', 'Average spacing: ***** Hz', 1, 'fontsize', 11, 'foregroundcolor', [0 0 0.25]);
P = GUIpanel('Zwuis', 'zwuis');
P = add(P, LowFreq, 'below', [20 0]);
P = add(P, HighFreq, nextto(LowFreq), [20 0]);
P = add(P, Ncomp, below(LowFreq));
P = add(P, BaseFreq, nextto(Ncomp), [10 0]);
P = add(P, FreqTol, alignedwith(Ncomp), [-5 0]);
P = add(P, Nit, nextto(FreqTol), [10 0]);
P = add(P, ZWseed, nextto(Nit), [5 0]);
P = add(P, ZwuisMess, below(FreqTol), [70 0]);

function P = local_profilePanel;
ProfileLowFreq = ParamQuery('ProfileLowFreq', 'low:', '15000.5', 'Hz', ...
    'rreal/nonnegative', 'Lowest frequency for profile sweep. Specify 0 Hz to always include lowest component', 1);
ProfileHighFreq = ParamQuery('ProfileHighFreq', 'high:', '15000.5', 'Hz', ...
    'rreal/positive', 'Highest frequency for profile sweep. Specify 100000 Hz to always include highest component', 1);
ProfileSPLjump = ParamQuery('ProfileSPLjump', 'SPL change:', '-15.5', 'dB', ...
    'rreal', 'Change in SPL of selected components.', 1);
ProfileType = ParamQuery('ProfileType', 'type:', '', {'single cmp' 'low end' 'high end' '-'}, ...
    '', 'How to realize profile. "-" = no profiling. Single cmp = SPL of one component is changed. Low/high: All components below/above swept freq are changed.', 1);
P = GUIpanel('Profile', 'profile');
P = add(P, ProfileLowFreq, 'below', [20 0]);
P = add(P, ProfileHighFreq, nextto(ProfileLowFreq), [20 0]);
P = add(P, ProfileSPLjump, below(ProfileLowFreq), [-20 0]);
P = add(P, ProfileType, nextto(ProfileSPLjump), [20 0]);

function P = local_presPanel;
TotDur = ParamQuery('TotDur', 'dur:', '150 ', 's', ...
    'rreal/positive', 'total duration per profile condition. Will be rounded to nearest integer # base cycles.',1);
Baseline = ParamQuery('Baseline', 'baseline:', '12000 12000 ', 'ms', ...
    'rreal/nonnegative', 'Duration of pre- and poststimulus baseline recording. Pairs are interpreted as [pre post].',2);
ISgap = ParamQuery('ISgap', 'gap:', '1500', 'ms', ...
    'rreal/positive', 'Gap betweenend of each stimulus presentation and the next one.',1);
RampDur = ParamQuery('RampDur', 'ramps:', '1500', 'ms', 'rreal/positive', 'Duration of ramps at on- & offset.',1);
Order = ParamQuery('Order', 'Order:', '', {'Forward' 'Reverse' 'Random'}, ...
    '', ['Play order of stimulus conditions. Forward means from Low to High value.' char(10) ...
    'Reverse means from High to Low. Random means conditions randomized.'],1);
RSeed = ParamQuery('RSeed', 'Rand Seed:', '844596300', '', ...
    'rseed', 'Random seed used for presentation order. Specify NaN to refresh seed upon each realization.',1);
P = GUIpanel('Pres', 'presentation');
P = add(P, TotDur,'below',[25 0]);
P = add(P,Baseline,nextto(TotDur), [0 0]);
P = add(P, ISgap,below(TotDur), [10 -5]);
P = add(P, RampDur,nextto(ISgap), [10 0]);
P = add(P,Order, alignedwith(ISgap), [0 -5]);
P = add(P,RSeed, nextto(Order), [0 0]);



