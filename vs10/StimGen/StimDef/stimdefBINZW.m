function Params = stimdefBINZW(EXP);
% stimdefBINZW - definition of stimulus and GUI for BINZW stimulus paradigm
%    P=stimdefBINZW(EXP) returns the definition for the BINZW (binzwuis)
%    stimulus paradigm. The definition P is a GUIpiece that can be rendered
%    by GUIpiece/draw. Stimulus definition like stimmdefZW are usually
%    called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimBINZW.

PairStr = ' Pairs of numbers are interpreted as [left right].';
ClickStr = ' Click button to select ';

ZW = local_zwuisPanel;
ITD = ITDstepper('ITD', EXP, '', 0); % last arg: specify different ITD types or not
SPL = SPLpanel('-', EXP, '' , 'Component');
Dur = DurPanel('Durations', EXP, '', 'basicsonly');
Pres = PresentationPanel;
%====================
Params=GUIpiece('Params'); % upper half of GUI: parameters
Params = add(Params, ZW, 'below', [20 0]);
Params = add(Params, ITD, below(ZW));
Params = add(Params, SPL, nextto(ZW), [20 0]);
Params = add(Params, Dur, nextto(ITD) ,[3 5]);
Params = add(Params, Pres, below(SPL) ,[200 2]);
Params = add(Params, PlayTime(), below(Dur) ,[0 5]);

%=================================
function P = local_zwuisPanel;
P = GUIpanel('ToneComplex', 'tones');
LowFreq = ParamQuery('LowFreq', 'low:', '15000.5', 'Hz', ...
    'rreal/positive', 'Approximate lowest frequency of tone complex.', 1);
HighFreq = ParamQuery('HighFreq', 'high:', '15000.5', 'Hz', ...
    'rreal/positive', 'Approximate highest frequency of tone complex.', 1);
Ncomp = ParamQuery('Ncomp', '#comp', '123', '', ...
    'posint', 'Number of tones in complex.' ,1);
FreqTol = ParamQuery('FreqTol', 'maxdev:', '120', '%', ...
    'rreal/positive', ['Maximum deviation of component frequencies from regular spacing, ' char(10), ...
    'expressed in % of average component spacing.'], 1);
ZWseed = ParamQuery('ZWseed', 'seed:', '12345678', '', ...
    'rseed', 'Random seed used to realize frequencies. Nan means refresh every time, EXP or EXP3, etc, means one seed for each experiment.' ,1);
ZwuisMess = messenger('ZwuisMess', 'Average spacing: ***** Hz', 1, 'fontsize', 11, 'foregroundcolor', [0 0 0.25]);
P = add(P, LowFreq, 'below', [20 0]);
P = add(P, HighFreq, nextto(LowFreq), [20 0]);
P = add(P, Ncomp, below(LowFreq));
P = add(P, FreqTol, alignedwith(Ncomp), [-5 0]);
P = add(P, ZWseed, nextto(FreqTol), [5 0]);
P = add(P, ZwuisMess, below(FreqTol), [70 0]);




