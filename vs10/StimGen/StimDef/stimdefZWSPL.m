function Params = stimdefZWSPL(EXP);
% stimdefZWSPL - definition of stimulus and GUI for ZWSPL stimulus paradigm
%    P=stimdefZWSPL(EXP) returns the definition for the ZWSPL (Zwuis SPL)
%    stimulus paradigm. The definition P is a GUIpiece that can be rendered
%    by GUIpiece/draw. Stimulus definition like stimmdefZW are usually
%    called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimZWSPL.

PairStr = ' Pairs of numbers are interpreted as [left right].';
ClickStr = ' Click button to select ';

ZW = local_zwuisPanel;
Dur = DurPanel('Durations', EXP, '', 'nophase', {'fine'});
SPL = SPLstepper('SPL', EXP); 
Pres = PresentationPanel;
%====================
Params=GUIpiece('Params'); % upper half of GUI: parameters
Params = add(Params, ZW, 'below', [20 0]);
Params = add(Params, Dur, below(ZW) ,[20 5]);
Params = add(Params, SPL, nextto(ZW), [10 0]);
Params = add(Params, Pres, below(SPL) ,[150 2]);
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




